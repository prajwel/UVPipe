# Use CentOS 6 as the base image
FROM centos:6.8

# Configure the vault repositories for CentOS 6
RUN sed -i 's|mirrorlist=|#mirrorlist=|g' /etc/yum.repos.d/CentOS-Base.repo && \
    sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-Base.repo

# Install necessary libraries
RUN yum install -y gcc gcc-c++ make cmake sqlite-devel csh

ENV AS1="/app"
ENV PFILES="$AS1/uvit/paramfiles:$PFILES"
ENV GLOG_log_dir="$AS1/uvit/log"
ENV LD_LIBRARY_PATH="$AS1/uvit/lib:$LD_LIBRARY_PATH"
ENV PATH="$AS1/uvit/bin:$PATH"

# Set the working directory in the container
WORKDIR /app

# Copy the entire codebase into the container
COPY . /app

# Run InstallLibs script
RUN chmod +x ./InstallLibs && ./InstallLibs

# Change to the uvit directory and build the project
WORKDIR /app/uvit
RUN make

WORKDIR /app/

CMD ["./VIS_multipipes.sh"]
