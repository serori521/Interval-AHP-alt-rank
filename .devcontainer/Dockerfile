FROM jupyter/datascience-notebook:julia-1.9.3

# Change this environment variable depending on the number of processers.
# See: https://docs.julialang.org/en/v1/manual/multi-threading/
ENV JULIA_NUM_THREADS=4

USER root

RUN groupadd -f jovyan && \
    usermod -a -G jovyan jovyan

RUN mkdir -p /home/jovyan/.cursor-server/data/User/globalStorage/julialang.language-julia/lsdepot/v1 && \
    chown -R jovyan:jovyan /home/jovyan/.cursor-server && \
    chmod -R u+rw /home/jovyan/.cursor-server

COPY ./init.jl /tmp/init.jl

USER jovyan

RUN julia /tmp/init.jl
