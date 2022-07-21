g++ -fopenmp \
-I /opt/local/include -L /opt/X11/lib -lX11 \
main.cpp -o main \
&& mkdir -p data \
&& mkdir -p data/con \
&& mkdir -p data/phi \
&& mkdir -p figures \
&& mkdir -p figures/con_xy \
&& mkdir -p figures/con_yz \
&& rm -f \
data/con/*.vtk \
data/con/*.csv \
data/phi/*.csv \
data/interface/*.csv \
data/fraction/*.csv \
figures/con_xy/*.png \
figures/con_yz/*.png \
figures/con/*.png \
figures/phi/*.png \
&& ./main \
&& rm main \
&& python plot1d.py
