## Authenticate with Earthdata Login

The first step to use NASA Earthdata is to create an account with Earthdata Login, please follow the instructions at [NASA EDL](https://urs.earthdata.nasa.gov/)

Once registered, earthaccess can use environment variables, a `.netrc` file or interactive input from a user to login with NASA EDL.

If a strategy is not especified, env vars will be used first, then netrc and finally user's input.

```py
import earthaccess

auth = earthaccess.login()
```

If you have a .netrc file with your Earthdata Login credentials

```py
auth = earthaccess.login(strategy="netrc")
```

If your Earthdata Login credentials are set as environment variables: EARTHDATA_USERNAME, EARTHDATA_PASSWORD

```py
auth = earthaccess.login(strategy="environment")
```

If you wish to enter your Earthdata Login credentials when prompted with optional persistence to .netrc

```py
auth = earthaccess.login(strategy="interactive", persist=True)
```
