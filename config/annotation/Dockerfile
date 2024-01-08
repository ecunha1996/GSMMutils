# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Update the package lists
RUN apt-get update && apt-get install -y wget

# Install OpenJDK 11
RUN apt-get install -y openjdk-11-jdk

# Install Python 3.7 and pip
RUN apt-get install -y python3.7 python3-pip

# Verify installations
RUN java -version
RUN python3 --version

# Get interproscan

RUN cd /home

#RUN wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.65-97.0/interproscan-5.65-97.0-64-bit.tar.gz
#RUN wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.65-97.0/interproscan-5.65-97.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
#RUN md5sum -c interproscan-5.65-97.0-64-bit.tar.gz.md5

#RUN tar -pxvzf interproscan-5.65-97.0-*-bit.tar.gz

#RUN python3 setup.py -f interproscan.properties