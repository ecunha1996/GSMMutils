
FROM ubuntu:22.04

ADD ./ /home/

RUN ls /home

RUN sh /home/cplex_studio2210.linux_x86_64.bin -f /home/cplex_response.txt

RUN apt-get update && apt-get install -y python3 python3-pip

RUN pip3 install --upgrade pip && pip3 install ray mewpy openpyxl

RUN cd /opt/ibm/ILOG/CPLEX_Studio221/cplex/python/3.10/x86-64_linux && python3 setup.py install && cd /home/