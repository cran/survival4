LIBNAME=survival4
LD=ld

OBJS= agexact.o      agsurv1.o      coxfit2.o      dmatrix.o      rnewton.o\
      agfit2.o       agsurv2.o      coxfit_null.o  doloop.o       \
      agfit_null.o   agsurv3.o      coxhaz2.o      survdiff2.o\
      aghaz2.o       chinv2.o       coxmart.o      pyears1.o      survexp2.o\
      agmart.o       cholesky2.o    coxres12.o     pyears2.o      survfit2.o\
      agres12.o      chsolve2.o     coxscho.o      pyears3.o      survindex2.o\
      agscore.o      coxdetail.o    coxscore.o     pystep.o       survreg.o   


survival4.so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so

realclean:
	@rm -f Makefile *.o *.so
