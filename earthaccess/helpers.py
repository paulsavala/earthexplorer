from ipyleaflet import Map, GeomanDrawControl, GeoJSON
from shapely.geometry import Polygon
from ipywidgets import Output
from IPython.display import display

def get_spatial_boundary():
    m = Map(center=(42, -90), zoom=3)

    draw_control = GeomanDrawControl()

    draw_control.marker = {
        "icon": {
            "iconUrl": "https://leafletjs.com/examples/custom-icons/leaf-green.png",
            "iconSize": [38, 95],
            "iconAnchor": [22, 94],
            "popupAnchor": [-3, -76],
            "shadowUrl": "https://leafletjs.com/examples/custom-icons/leaf-shadow.png",
            "shadowSize": [50, 64],
            "shadowAnchor": [4, 62],
        },
    }

    draw_control.polygon = {
        "pathOptions": {
            "fillColor": "blue",
            "color": "blue",
            "fillOpacity": 0.25,
        }
    }

    draw_control.polyline = {}
    draw_control.circlemarker = {}
    draw_control.cut = False
    draw_control.remove = False
    draw_control.rotate = False
    draw_control.edit = False
    draw_control.remove = True

    m.add(draw_control)

    # Create an output widget to display results
    out = Output()

    # Define callback function for draw events
    def handle_draw(self, action, geo_json):
        with out:
            if action == 'remove':
                out.clear_output()
            if action in ['create', 'drag']:
                out.clear_output()
                geo_json = geo_json[0]
                if geo_json['geometry']['type'] == 'Point':
                    lon, lat = geo_json['geometry']['coordinates']
                    print([lat, lon])
                    return (lat, lon)
                elif geo_json['geometry']['type'] == 'Polygon':  # Rectangle is a type of Polygon
                    coords = geo_json['geometry']['coordinates'][0]
                    print('[')
                    for lon, lat in coords:
                        print(f'  [{lat}, {lon}],')
                    print(']')
                    return coords

    # Attach the callback to the DrawControl
    draw_control.on_draw(handle_draw)

    # Display the map and output widget
    display(m, out)
    
    
def display_spatial_boundary(coords):
    if isinstance(coords[0], list):
        polygon = Polygon(coords)
        
        # Check if the polygon is valid
        if not polygon.is_valid:
            polygon = polygon.buffer(0)

        # Ensure counterclockwise orientation
        if not polygon.exterior.is_ccw:
            polygon = Polygon(list(polygon.exterior.coords)[::-1])

        coords = list(polygon.exterior.coords)

        # Ensure the polygon is closed
        if coords[0] != coords[-1]:
            coords.append(coords[0])
            
    if coords != []:
        if isinstance(coords[0], float):
            polygon_geojson = {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "Point",
                            "coordinates": [coords[-1], coords[0]]  # GeoJSON uses [lon, lat] format
                        },
                    }
                ]
            }
        else:
            polygon_geojson = {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [[(lon, lat) for lat, lon in coords]]  # GeoJSON uses [lon, lat] format
                        },
                        "properties": {
                            "index": 1 # An index is required, otherwise it won't display
                        }
                    }
                ]
            }
            
        m = Map(center=[42, -90], zoom=3)
        
        # Create a GeoJSON layer with the polygon
        polygon_layer = GeoJSON(data=polygon_geojson)
        m.add_layer(polygon_layer)
        
        return m
