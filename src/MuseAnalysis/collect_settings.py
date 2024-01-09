import glob
import json
import pandas as pd
import argparse


def main():
    """
    Collect all the setting files in the subfolders and put them into a single json and csv file
    """
    
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("-p", "--prefix")
    parser.add_argument("-o","--outputname")
    
    args = parser.parse_args()
    if hasattr(args, "help"):
        if args.help:
            parser.print_help()
    if hasattr(args, "prefix"):
        prefix = args. prefix
    else:
        prefix = ""
    if hasattr(args, "outputname"):
        out = args.outputname
    else:
        out = "log"
    
    shots = glob.glob(prefix+"*[0-9]")
    shots.sort()

    logs: dict[str, dict] = {}

    for shot in shots:
        with open(shot + "/settings.txt") as f:
            lines = f.readlines()

        logs[shot] = {}
        for line in lines:
            head, body = line.split(":")
            head = head.strip()
            body = body.strip()
            if head in logs[shot]:
                head = head + "_1"
            logs[shot][head] = body

    with open(out+".json", "w") as f:
        json.dump(logs, f)

    df = pd.DataFrame(data=logs).T
    df.to_csv(out+".csv", sep="\t")

    return logs


if __name__ == "__main__":
    logs = main()
    print(logs)
