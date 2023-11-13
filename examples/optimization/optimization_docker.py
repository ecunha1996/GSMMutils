from GSMMutils.optimization.optimization_docker import OptimizationDocker


def main():
    optimization = OptimizationDocker(data_directory="/home/ecunha/tmp/pycharm_project_27/data",
                                      src_directory="/home/ecunha/tmp/pycharm_project_27/src",
                                      examples_directory="/home/ecunha/tmp/pycharm_project_27/examples",
                                      config_directory="/home/ecunha/tmp/pycharm_project_27/config",
                                      utilities_directory="/home/ecunha/tmp/pycharm_project_27/utilities")
    # optimization.build()
    optimization.run()


if __name__ == '__main__':
    main()
