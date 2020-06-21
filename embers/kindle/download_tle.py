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
