from oorinano.simulator.msdft import draw_QFL
import argparse
import os, sys

def main():
    parser = argparse.ArgumentParser(description="run qfl visualization with data")
    parser.add_argument('-d', '--dir', default='QFL_data',  help='directory with wavefunction data')
    parser.add_argument('-u', '--usage', action='store_true', help='usage for run_qt')
    args = parser.parse_args()

    if args.usage:
        print(f"Usage::\
            python visual_QFL.py -d {args.dir}")
        sys.exit(1)
    print(f"QFL data path {args.dir}")
    draw_QFL(args.dir)

if __name__ == "__main__":
    main()
