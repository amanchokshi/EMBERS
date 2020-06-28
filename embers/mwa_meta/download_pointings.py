import json, wget, time
from pathlib import Path
from astropy.time import Time


def wget_meta(start, stop, num_pages, out_dir):
    start_gps = int(Time(start, format="isot").gps)
    stop_gps = int(Time(stop, format="isot").gps)
    for npg in range(num_pages):
        time.sleep(28)
        cerberus_url = f"http://ws.mwatelescope.org/metadata/find?mintime={start_gps}&maxtime={stop_gps}&extended=1&page={npg+1}&pretty=1"
        print(f"\nDownloading page {npg+1} of metadata")
        wget.download(cerberus_url, f"{out_dir}/pointing_{npg+1:03d}.json")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="""
            Download pointing metadata
            """
    )

    parser.add_argument(
        "--start",
        metavar="\b",
        help="Start date in isot format. Ex: 2019-10-01T00:00:00",
    )
    parser.add_argument(
        "--stop",
        metavar="\b",
        help="Stop date in isot format. Ex: 2019-10-10T00:00:00",
    )
    parser.add_argument(
        "--num_pgs",
        type=int,
        metavar="\b",
        help="Total number of pages between start and stop, determined from http://ws.mwatelescope.org/metadata/find",
    )
    parser.add_argument(
        "--out_dir", metavar="\b", help="Output directory for json metadata",
    )

    args = parser.parse_args()
    start = args.start
    stop = args.stop
    num_pgs = args.num_pgs
    out_dir = Path(args.out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    wget_meta(start, stop, num_pgs, out_dir)
