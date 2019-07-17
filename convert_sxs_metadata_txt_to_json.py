import argparse
import sxs

p = argparse.ArgumentParser(description="Convert a metadata.txt file to JSON")
p.add_argument("--txt", help="Path to input metadata.txt", required=True)
p.add_argument("--output", help="Path to write JSON file", required=True)
args = p.parse_args()

sxs.metadata.Metadata.from_file(args.txt).to_json_file(args.output)
