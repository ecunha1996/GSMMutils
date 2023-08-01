from GSMMutils.omics.troppo_docker import TroppoDocker


def main():
    dfba = TroppoDocker(data_directory="/home/ecunha/tmp/pycharm_project_27/data", src_directory="/home/ecunha/tmp/pycharm_project_27/src", examples_directory="/home/ecunha/tmp/pycharm_project_27/examples",
                        config_directory="/home/ecunha/tmp/pycharm_project_27/config")
    dfba.run()


if __name__ == '__main__':
    main()
