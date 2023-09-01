import sys
sys.path.append('../../src')
from GSMMutils.dynamic.dfba_docker import DFBA


def main():
    dfba = DFBA(data_directory="/home/ecunha/tmp/pycharm_project_27/data", src_directory="/home/ecunha/tmp/pycharm_project_27/src")
    dfba.run()


if __name__ == '__main__':
    main()
