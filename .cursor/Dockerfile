FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# System image only. Cursor checks out the commit; do not COPY the project.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    gnupg \
    nodejs \
    openjdk-21-jdk-headless \
    unzip \
    zip \
    && rm -rf /var/lib/apt/lists/*

# sbt launcher; the project pin is project/build.properties (1.12.14).
RUN echo "deb https://repo.scala-sbt.org/scalasbt/debian all main" \
      > /etc/apt/sources.list.d/sbt.list \
    && curl -fsSL "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x2EE0EA64E40A89B84B2DF73499E82A75642AC823" \
      | gpg --dearmor -o /etc/apt/trusted.gpg.d/sbt.gpg \
    && apt-get update \
    && apt-get install -y --no-install-recommends sbt \
    && rm -rf /var/lib/apt/lists/*
