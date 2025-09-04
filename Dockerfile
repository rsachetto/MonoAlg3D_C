
FROM ubuntu:plucky

RUN apt-get update && apt-get install -y \
    build-essential \
    make \
    libglu1-mesa-dev \
    libx11-dev \
    libopengl-dev \
    libglx-dev \
    libxcursor-dev \
    libxi-dev \
    libxrandr-dev \
    libxinerama-dev \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/monoalg
WORKDIR /opt/monoalg
RUN make clean
RUN make simulator || true
RUN ./build.sh gui || true
CMD ["bash"]
