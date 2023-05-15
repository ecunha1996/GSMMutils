
import sys
sys.path.append('../../src')

from ExpGSMM.dynamic.dfba_docker import dFBA


def main():
    dfba = dFBA(data_directory="/home/ecunha/tmp/pycharm_project_27/data", src_directory="/home/ecunha/tmp/pycharm_project_27/src")
    dfba.run()



if __name__ == '__main__':
    main()







