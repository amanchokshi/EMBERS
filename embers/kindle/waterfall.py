import argparse
from embers.condition_data.rf_data import single_waterfall



parser = argparse.ArgumentParser(
    description="""
    Saved a single waterfall plot
    """
)

parser.add_argument("--rf_file", metavar="\b", help="Path to raw rf data file")
parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/condition_data",
    help="Dir where colormap sample plot is saved. Default=./embers_out/condition_data",
)

args = parser.parse_args()
rf_file = args.rf_file
out_dir = args.out_dir

single_waterfall(rf_file, out_dir)


