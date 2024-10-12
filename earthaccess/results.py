import json
import uuid
import geopandas as gpd
import pandas as pd
import numpy as np
from typing import Any, Dict, List, Tuple, Optional, Union, cast, Type
from ipyleaflet import Map, GeoJSON, DrawControl
from itables import show
from datetime import datetime, timezone
from matplotlib import pyplot as plt
from shapely.geometry import Polygon, Point

from .formatters import _repr_granule_html
from typing import Any, Dict, List, Optional, Union
from .services import DataServices
from .auth import Auth
from .utils._results import convert_to_mb

import earthaccess


class CustomDict(dict):
    _basic_umm_fields_: List = []
    _basic_meta_fields_: List = []

    def __init__(
        self,
        collection: Dict[str, Any],
        fields: Optional[List[str]] = None,
        cloud_hosted: bool = False,
    ):
        super().__init__(collection)
        self.cloud_hosted = cloud_hosted
        self.uuid = str(uuid.uuid4())

        self.render_dict: Any
        if fields is None:
            self.render_dict = self
        elif fields[0] == "basic":
            self.render_dict = self._filter_fields_(self._basic_umm_fields_)
        else:
            self.render_dict = self._filter_fields_(fields)

    def _filter_fields_(self, fields: List[str]) -> Dict[str, Any]:
        filtered_dict = {
            "umm": dict(
                (field, self["umm"][field]) for field in fields if field in self["umm"]
            )
        }
        basic_dict = {
            "meta": dict(
                (field, self["meta"][field])
                for field in self._basic_meta_fields_
                if field in self["meta"]
            )
        }
        basic_dict.update(filtered_dict)
        return basic_dict

    def _filter_related_links(self, filter: str) -> List[str]:
        """Filter RelatedUrls from the UMM fields on CMR."""
        matched_links: List = []
        if "RelatedUrls" in self["umm"]:
            for link in self["umm"]["RelatedUrls"]:
                if link["Type"] == filter:
                    matched_links.append(link["URL"])
        return matched_links




class DataCollection(CustomDict):
    """Dictionary-like object to represent a data collection from CMR."""

    _basic_meta_fields_ = [
        "concept-id",
        "granule-count",
        "provider-id",
    ]

    _basic_umm_fields_ = [
        "ShortName",
        "Abstract",
        "SpatialExtent",
        "TemporalExtents",
        "DataCenters",
        "RelatedUrls",
        "ArchiveAndDistributionInformation",
        "DirectDistributionInformation",
    ]

    def summary(self) -> Dict[str, Any]:
        """Summary containing short_name, concept-id, file-type, and cloud-info (if cloud-hosted).

        Returns:
            A summary of the collection metadata.
        """
        # we can print only the concept-id

        summary_dict: Dict[str, Any]
        summary_dict = {
            "short-name": self.get_umm("ShortName"),
            "entity-title": self.get_umm("EntryTitle"),
            "concept-id": self.concept_id(),
            "version": self.version(),
            "file-type": self.data_formats(),
            "size": self.size(),
            "num-granules": self.get_meta("granule-count"),
            "get-data": self.landing_page(),
        }
        if "Region" in self.s3_bucket():
            summary_dict["cloud-info"] = self.s3_bucket()
        return summary_dict

    def get_umm(self, umm_field: str) -> Union[str, Dict[str, Any]]:
        """Placeholder.

        Parameters:
            umm_field: Valid UMM item, i.e. `TemporalExtent`.

        Returns:
            The value of a given field inside the UMM (Unified Metadata Model).
        """
        if umm_field in self["umm"]:
            return self["umm"][umm_field]
        return ""
    
    def get_meta(self, meta_field: str) -> Union[str, Dict[str, Any]]:
        """
        Parameters:
            meta_field: Valid meta item, i.e. `granule-count`

        Returns:
            The value of a given field inside the meta (collection metadata).
        """
        if meta_field in self["meta"]:
            return self["meta"][meta_field]
        return ""
    
    def size(self) -> Optional[float]:
        """
        Returns:
            The total size for the granule in MB.
        """
        try:
            data_collection = self["umm"]
            data_collection_info = data_collection["ArchiveAndDistributionInformation"]
            total_size = sum(
                [
                    convert_to_mb(float(s["TotalCollectionFileSize"]), s["TotalCollectionFileSizeUnit"])
                    for s in data_collection_info["FileDistributionInformation"]
                    if "FileDistributionInformation" in data_collection_info
                ]
            )
        except Exception as e:
            print(e)
            try:
                data_collection = self["umm"]
                data_collection_info = data_collection["ArchiveAndDistributionInformation"]
                total_size = sum(
                    [
                        float(s["SizeInBytes"])
                        for s in data_collection_info["FileDistributionInformation"]
                        if "FileDistributionInformation" in data_collection_info
                    ]
                ) / (1024 * 1024)
            except Exception as e:
                total_size = 0.
        return total_size
    
    def data_formats(self) -> str:
        """
        Returns:
            The collection data type, i.e. HDF5, CSV etc., if available, as a comma-separated string.
        """
        if "ArchiveAndDistributionInformation" in self["umm"]:
            file_format_set = set()
            for file in self["umm"]["ArchiveAndDistributionInformation"]["FileDistributionInformation"]:
                file_format_set.add(str(file['Format']))
                    
            return ", ".join(file_format_set)
        return ""

    def concept_id(self) -> str:
        """
        Returns:
            A collection's `concept_id`.This id is the most relevant search field on granule queries.
        """
        return self["meta"]["concept-id"]

    def version(self) -> str:
        """

        Returns:
            The collection's version.
        """
        if "Version" in self["umm"]:
            return self["umm"]["Version"]
        return ""

    def abstract(self) -> str:
        """Placeholder.

        Returns:
            The abstract of a collection.
        """
        if "Abstract" in self["umm"]:
            return self["umm"]["Abstract"]
        return ""

    def landing_page(self) -> str:
        """Placeholder.

        Returns:
            The first landing page for the collection (can be many), if available.
        """
        links = self._filter_related_links("LANDING PAGE")
        if len(links) > 0:
            return links[0]
        return ""
    
    def get_horizontal_resolution(self) -> None:
        pass

    def get_data(self) -> List[str]:
        """
        Returns:
            The GET DATA links (usually a landing page link, a DAAC portal, or an FTP location).
        """
        links = self._filter_related_links("GET DATA")
        return links

    def s3_bucket(self) -> Dict[str, Any]:
        """Placeholder.

        Returns:
            The S3 bucket information if the collection has it (**cloud hosted collections only**).
        """
        if "DirectDistributionInformation" in self["umm"]:
            return self["umm"]["DirectDistributionInformation"]
        return {}

    def services(self) -> Dict[Any, List[Dict[str, Any]]]:
        """Return list of services available for this collection."""
        services = self.get("meta", {}).get("associations", {}).get("services", [])
        queries = (
            DataServices(auth=cast(Auth, earthaccess.__auth__)).parameters(concept_id=service)
            for service in services
        )

        return {service: query.get_all() for service, query in zip(services, queries)}

    def __repr__(self) -> str:
        return json.dumps(
            self.render_dict, sort_keys=False, indent=2, separators=(",", ": ")
        )


class DataCollectionList(List):
    """List-like object to represent a list of collections from CMR."""

    def __init__(self, granules):
        super().__init__(granules)

    # def __repr__(self) -> List:
    #     """
    #     Returns:
    #         The normal list of all collections (i.e. acts just like a list).
    #     """
    #     return super().__repr__()
    
    def _create_button_html(self, link):
        return f'<button onclick="window.open(\'{link}\', \'_blank\')">View Data</button>'
    
    def summary(self, show_links=False) -> Union[pd.DataFrame, None]:
        """
        Create a summary table of search results.
        Parameters
        ----------
        show_links : bool
            If True, the table will include a column with buttons to view the data.

        Returns
        -------
        pd.DataFrame
            Summary table of search results.
        """
        summary_list = [result.summary() for result in self]
        summary_df = pd.DataFrame(summary_list)

        # Change cloud hosting column to a bool
        summary_df['cloud-hosted'] = summary_df['cloud-info'].map(lambda x: 'Region' in x.keys())
        summary_df = summary_df.drop(columns='cloud-info')

        # Convert data links to buttons
        if show_links:
            buttons_html = [self._create_button_html(link) for link in summary_df['get-data']]

            # Insert the buttons HTML as a new column in the DataFrame
            summary_df['view-data'] = buttons_html
            summary_df = summary_df.drop(columns='get-data')

        # Reorder columns to put most useful columns first
        col_order = ['concept-id', 'entity-title', 'short-name', 'file-type', 'size', 'num-granules', 'cloud-hosted']
        
        if show_links:
            col_order = ['view-data'] + col_order

        summary_df = summary_df[col_order]

        if show_links:
            show(summary_df)
        else:
            return summary_df


class DataGranuleList(List):
    """List-like object to represent a list of granules from CMR."""

    def __init__(self, granules):
        super().__init__(granules)
        
    def __len__(self) -> int:
        return super().__len__()

    # def __repr__(self) -> List:
    #     """
    #     Returns:
    #         The normal list of all granules (i.e. acts just like a list).
    #     """
    #     return super().__repr__()
    
    def summary(self) -> pd.DataFrame:
        """
        Create a summary table of granules.

        Parameters
        ----------
        results : List
            List of search results from `earthaccess.search_datasets` or `earthaccess.search_data`.

        Returns
        -------
        pd.DataFrame
            Summary table of search results.
        """
        summary_list = [result.summary() for result in self]
        summary_df = pd.DataFrame(summary_list)

        # Reorder columns to put most useful columns first
        col_order = ['concept-id', 'file-size', 'num-files', 'start-datetime', 'end-datetime']
        summary_df = summary_df[col_order]
        summary_df['start-datetime'] = pd.to_datetime(summary_df['start-datetime'])
        summary_df['end-datetime'] = pd.to_datetime(summary_df['end-datetime'])

        return summary_df

    def _datetime_range_within(self, dt_range: Tuple[datetime], start: Optional[datetime] = None, finish: Optional[datetime] = None):
        # Extract the start and end times
        dt_start, dt_finish = dt_range[0], dt_range[-1]
        
        # Check if the dt_range is within the start and finish
        # Default to True if no start or finish is supplied (no range to check against)
        within_range = True
        if start is not None:
            within_range = within_range and dt_start >= start
        if finish is not None:
            within_range = within_range and dt_finish <= finish
        
        return within_range

    def filter_temporal_extent(self, start_dt: Optional[Union[str, datetime]] = None, finish_dt: Optional[Union[str, datetime]] = None) -> 'DataGranuleList':
        start_dt_utc = None
        finish_dt_utc = None

        if start_dt is not None:
            if isinstance(start_dt, str):
                start_dt_utc = datetime.fromisoformat(start_dt).replace(tzinfo=timezone.utc)
            else:
                start_dt_utc = start_dt.replace(tzinfo=timezone.utc)
        if finish_dt is not None:
            if isinstance(finish_dt, str):
                finish_dt_utc = datetime.fromisoformat(finish_dt).replace(tzinfo=timezone.utc)
            else:
                finish_dt_utc = finish_dt.replace(tzinfo=timezone.utc)

        filtered_granules = [granule for granule in self 
                            if self._datetime_range_within(dt_range=granule.temporal_extent(), 
                                                            start=start_dt_utc, 
                                                            finish=finish_dt_utc)]

        return DataGranuleList(filtered_granules)
    
    def plot_temporal_extent(self):
        """
        Plots a histogram of all temporal extents in the DataGranuleList.
        """
        start_dates = [granule.temporal_extent()[0] for granule in self]
        end_dates = [granule.temporal_extent()[1] for granule in self]
        
        # Combine start and end dates for histogram
        all_dates = start_dates + end_dates

        ax = plt.figure().gca()
        plt.hist(all_dates, bins=50, alpha=0.7, color='blue', edgecolor='black')
        plt.title('Histogram of Temporal Extents')
        plt.xlabel('Date')
        plt.ylabel('Frequency')
        ax.yaxis.get_major_locator().set_params(integer=True) # Only show integer ticks
        plt.show()

    def plot_spatial_extents(self):
        granule_extents = []
        granules = []
        for granule in self:
            if granule.get_spatial_extent() not in granule_extents:
                granule_extents.append(granule.get_spatial_extent())
                granules.append(granule)
        
        m = Map(center=[42, -90], zoom=3)
        for granule in granules:
            spatial_extent = granule.get_spatial_extent()
            concept_id = granule.get_meta('concept-id')
            layer = GeoJSON(data=granule._create_polygon_geojson(spatial_extent, concept_id))
            m.add_layer(layer)
        
        return m

    def filter_spatial_extents(self, coords: Union[List, List[List[tuple[float, float]]]], method: str) -> 'DataGranuleList':
        """
        Filter granules based on spatial extent overlap with the drawn rectangle.
        
        Parameters
        ----------
        coords : List of [longitude, latitude] coordinates defining the rectangle
        
        Returns
        -------
        DataGranuleList : A filtered list of granules with spatial overlap
        """
        method = method.lower().strip()
        if not method in ['intersect', 'contains', 'within']:
            raise ValueError('Invalid method. Please choose from "intersect", "contains", or "within".')

        is_polygon = False
        if isinstance(coords[0], float):
            coords_point = Point(coords)
        else:
            coords_polygon = Polygon(coords)
            is_polygon = True
            
        if not is_polygon:
            if method != '' and method != 'within':
                raise ValueError('Point most use the "within" method. Either provide a polygon or change the method.')
            else:
                method = 'within'
        
        filtered_granules = []
        for granule in self:
            # Get the spatial extent of the granule (assuming it's a polygon)
            granule_coords = granule.get_spatial_extent()
            for granule_coord in granule_coords:
                granule_polygon = Polygon(granule_coord)
                
                if method == 'intersect' and granule_polygon.intersects(coords_polygon):
                    filtered_granules.append(granule)
                if method == 'contains' and coords_polygon.contains(granule_polygon):
                    filtered_granules.append(granule)
                if method == 'within':
                    if is_polygon and granule_polygon.contains(coords_polygon):
                        filtered_granules.append(granule)
                    elif not is_polygon and granule_polygon.contains(coords_point):
                        filtered_granules.append(granule)        
        
        return DataGranuleList(filtered_granules)
    

class DataGranule(CustomDict):
    """Dictionary-like object to represent a granule from CMR."""

    _basic_meta_fields_ = [
        "concept-id",
        "provider-id",
    ]

    _basic_umm_fields_ = [
        "GranuleUR",
        "SpatialExtent",
        "TemporalExtent",
        "RelatedUrls",
        "DataGranule",
    ]

    def __init__(
        self,
        collection: Dict[str, Any],
        fields: Optional[List[str]] = None,
        cloud_hosted: bool = False,
    ):
        super().__init__(collection)
        self.cloud_hosted = cloud_hosted
        # TODO: maybe add area, start date and all that as an instance value
        self["size"] = self.size()
        self.uuid = str(uuid.uuid4())
        self.render_dict: Any
        if fields is None:
            self.render_dict = self
        elif fields[0] == "basic":
            self.render_dict = self._filter_fields_(self._basic_umm_fields_)
        else:
            self.render_dict = self._filter_fields_(fields)

    def __repr__(self) -> str:
        """Placeholder.

        Returns:
            A basic representation of a data granule.
        """
        data_links = [link for link in self.data_links()]
        rep_str = f"""
        Collection: {self['umm'].get('CollectionReference', "N/A")}
        Spatial coverage: {self['umm'].get('SpatialExtent', "N/A")}
        Temporal coverage: {self['umm'].get('TemporalExtent', "NA")}
        Size(MB): {self.size():.2f}
        Data: {data_links}\n\n
        """.strip().replace("  ", "")
        return rep_str

    def _repr_html_(self) -> str:
        """Placeholder.

        Returns:
            A rich representation for a data granule if we are in a Jupyter notebook.
        """
        granule_html_repr = _repr_granule_html(self)
        return granule_html_repr
    
    def summary(self) -> Dict[str, Any]:
        """Summary containing short_name, concept-id, file-type, and cloud-info (if cloud-hosted).

        Returns:
            A summary of the collection metadata.
        """

        summary_dict: Dict[str, Any]
        summary_dict = {
            "concept-id": self.get_meta('concept-id'),
            "file-size": self.size(),
            "num-files": self.num_files(),
            "start-datetime": self.temporal_extent()[0],
            "end-datetime": self.temporal_extent()[1],
        }
        return summary_dict
    
    def get_umm(self, umm_field: str) -> Union[str, Dict[str, Any]]:
        """
        Parameters:
            umm_field: Valid UMM item, i.e. `TemporalExtent`

        Returns:
            The value of a given field inside the UMM (Unified Metadata Model).
        """
        if umm_field in self["umm"]:
            return self["umm"][umm_field]
        return ""
    
    def get_meta(self, meta_field: str) -> Union[str, Dict[str, Any]]:
        """
        Parameters:
            meta_field: Valid meta item, i.e. `granule-count`

        Returns:
            The value of a given field inside the meta (collection metadata).
        """
        if meta_field in self["meta"]:
            return self["meta"][meta_field]
        return ""

    def get_s3_credentials_endpoint(self) -> Optional[str]:
        for link in self["umm"]["RelatedUrls"]:
            if "/s3credentials" in link["URL"]:
                return link["URL"]
        return None

    def size(self) -> float:
        """
        Returns:
            The total size for the granule in MB.
        """
        total_size = 0.
        try:
            data_granule = self["umm"]["DataGranule"]
            if "ArchiveAndDistributionInformation" in data_granule:
                data_granule_info = data_granule["ArchiveAndDistributionInformation"]
                total_size = sum(
                    [
                        convert_to_mb(float(s["Size"]), s["SizeUnit"])
                        for s in data_granule_info
                        if "Size" in s.keys()
                    ]
                )
        except Exception:
            try:
                data_granule = self["umm"]["DataGranule"]
                total_size = sum(
                    [
                        convert_to_mb(float(s["SizeInBytes"]), 'BYTES')
                        for s in data_granule["ArchiveAndDistributionInformation"]
                        if "ArchiveAndDistributionInformation" in data_granule
                    ]
                ) / (1024 * 1024)
            except Exception:
                total_size = 0.
        return total_size
    
    def num_files(self) -> int:
        """
        Returns:
            The number of files inside this granule.
        """
        if "ArchiveAndDistributionInformation" in self["umm"]['DataGranule']:
            files = [x for x in 
                     self['umm']['DataGranule']['ArchiveAndDistributionInformation']
                     if 'Size' in x.keys()]
            num_files = len(files)
            return num_files
        return 0

    def _derive_s3_link(self, links: List[str]) -> List[str]:
        s3_links = []
        for link in links:
            if link.startswith("s3"):
                s3_links.append(link)
            elif link.startswith("https://") and (
                "cumulus" in link or "protected" in link
            ):
                s3_links.append(f's3://{links[0].split("nasa.gov/")[1]}')
        return s3_links

    def data_links(
        self, access: Optional[str] = None, in_region: bool = False
    ) -> List[str]:
        """Returns the data links from a granule.

        Parameters:
            access: direct or external.
                Direct means in-region access for cloud-hosted collections.
            in_region: True if we are running in us-west-2.
                It is meant for the store class.

        Returns:
            The data links for the requested access type.
        """
        https_links = self._filter_related_links("GET DATA")
        s3_links = self._filter_related_links("GET DATA VIA DIRECT ACCESS")
        if in_region:
            # we are in us-west-2
            if self.cloud_hosted and access in (None, "direct"):
                # this is a cloud collection, and we didn't specify the access type
                # default to S3 links
                if len(s3_links) == 0 and len(https_links) > 0:
                    # This is guessing the S3 links for some cloud collections that for
                    # some reason only offered HTTPS links
                    return self._derive_s3_link(https_links)
                else:
                    # we have the s3 links so we return those
                    return s3_links
            else:
                # Even though we are in us-west-2, the user wants the HTTPS links used in-region.
                # They are S3 signed links from TEA.
                # <https://github.com/asfadmin/thin-egress-app>
                return https_links
        else:
            # we are not in-region
            if access == "direct":
                # maybe the user wants to collect S3 links and use them later
                # from the cloud
                return s3_links
            else:
                # we are not in us-west-2, even cloud collections have HTTPS links
                return https_links

    def dataviz_links(self) -> List[str]:
        """
        Returns:
            The data visualization links, usually the browse images.
        """
        links = self._filter_related_links("GET RELATED VISUALIZATION")
        return links
    
    def _create_polygon_geojson(self, coords, granule_id) -> Dict[str, Any]:
        """
        Returns:
            A GeoJSON object with the spatial extent of the granule.
        """
        polygon_geojson = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[(lat, lon) for lat, lon in sublist]]  # GeoJSON uses [lon, lat] format
                    },
                    "properties": {
                        "index": granule_id
                    }
                } for sublist in coords
            ]
        }
        return polygon_geojson
    
    def get_spatial_extent(self, lat_first: bool = True) -> List[List[tuple[float, float]]]:
        """
        Returns:
            A list of tuples of floats showing the spatial extent for each granule.
        """
        try:
            spatial_extent_polygons = self['umm']['SpatialExtent']['HorizontalSpatialDomain']['Geometry']['GPolygons']
            coords = [
                        [
                            (point['Latitude'], point['Longitude']) if lat_first else (point['Longitude'], point['Latitude'])
                            for point in boundary['Boundary']['Points']
                        ]
                        for boundary in spatial_extent_polygons
                    ]
        except KeyError:
            coords = []
        return coords
    
    def plot_spatial_extent(self) -> Union[Map, None]:
        """
        Returns:
            An ipyleaflet map with the spatial extent of the granule.
        """
        # Create a starter map, centered on North America
        m = Map(center=[42, -90], zoom=3)
        
        # Note: GeoJSON expects longitude first, then latitude: https://geojson.org/geojson-spec.html#positions
        coords = self.get_spatial_extent(lat_first=False)
        
        if coords != []:
            polygon_geojson = self._create_polygon_geojson(coords, self.get_meta('concept-id'))
            
            # Create a GeoJSON layer with the polygon
            polygon_layer = GeoJSON(data=polygon_geojson)
            m.add_layer(polygon_layer)
        else:
            print('No spatial extent found for this granule.')
            m = None

        # Return the map
        return m
    
    def temporal_extent(self) -> tuple[datetime, datetime]:
        """
        Returns:
            The temporal extent of the granule, as a tuple with start and end times, respectively.
        """
        start_dt = self['umm']['TemporalExtent']['RangeDateTime']['BeginningDateTime']
        end_dt = self['umm']['TemporalExtent']['RangeDateTime']['EndingDateTime']

        return (datetime.fromisoformat(start_dt).replace(tzinfo=timezone.utc), 
                datetime.fromisoformat(end_dt).replace(tzinfo=timezone.utc))
