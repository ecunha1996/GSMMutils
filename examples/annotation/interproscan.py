from GSMMutils.annotation.interproscan import InterProScanDocker


def main():
    interproscan = InterProScanDocker(config_directory="/home/ecunha/tmp/pycharm_project_27/config")
    interproscan.build()






if __name__ == '__main__':
    main()