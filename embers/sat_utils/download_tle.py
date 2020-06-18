import os
import time
import argparse
from pathlib import Path
from sat_ids import norad_ids
import spacetrack.operators as op
from spacetrack import SpaceTrackClient

# space-tracks.org login and password, saved as environment variables

if os.environ.get("ST_USER") != None:
    st_user = os.environ.get("ST_USER")
    st_pass = os.environ.get("ST_PASS")

    st = SpaceTrackClient(identity=st_user, password=st_pass)

    parser = argparse.ArgumentParser(
        description="""Downloads TLE files for the Beam Experiment from Space-tracks.org."""
    )
    parser.add_argument(
        "--start_date", metavar="\b", help="Start date in format YYYY-MM-DD"
    )
    parser.add_argument(
        "--stop_date", metavar="\b", help="Stop date in format YYYY-MM-DD"
    )
    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./../../outputs/sat_ephemeris/TLE/",
        help="Path to output directory. Default=./../../outputs/sat_ephemeris/TLE/",
    )

    args = parser.parse_args()

    start_date = args.start_date
    stop_date = args.stop_date
    out_dir = args.out_dir

    # make a TLE directory
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    os.makedirs(os.path.dirname(out_dir), exist_ok=True)

    print("Starting TLE download")

    for sat_name, sat_id in norad_ids.items():
        # Sleep to limit downloads to 200 TLEs per hour
        time.sleep(20)
        data = st.tle(
            iter_lines=True,
            norad_cat_id=sat_id,
            orderby="epoch asc",
            epoch=f"{start_date}--{stop_date}",
            format="tle",
        )
        print(
            f"downloading tle for {sat_name} satellite [{sat_id}] from space-tracks.org"
        )
        with open(f"{out_dir}/{sat_id}.txt", "w") as fp:
            for line in data:
                fp.write(line + "\n")


else:
    print("ERROR!")
    print("Space-Track.org credentials do not exist!")
    print("Edit this file and enter your account details.")
