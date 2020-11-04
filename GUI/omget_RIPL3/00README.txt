
gfortran -c omget.f om_retrieve.f
gfortran -o omget.exe omget.o om_retrieve.o -static -static-libgfortran