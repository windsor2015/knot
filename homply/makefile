CFLAGS = -O2

.cc.o:
	gcc $(CFLAGS) -c $<

O = dllink.o poly.o knot.o bound.o order.o model.o control.o

homfly : $(O)
	ar rcs libhomfly.a $(O)

# DEPENDENCIES

dllink.o : standard.h dllink.h dllink.c

poly.o : standard.h poly.c

knot.o : standard.h dllink.h knot.h

bound.o : standard.h order.h control.h bound.h

order.o : standard.h dllink.h knot.h bound.h order.h

model.o : standard.h poly.h order.h bound.h control.h model.h

control.o : standard.h poly.h knot.h order.h bound.h model.h control.h

clean :
	rm *.o
