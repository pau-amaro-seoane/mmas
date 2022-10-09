NAME = meastar
NAME2= table2
NAME3= fig25

FC= egfortran
FFILES = driver.f makemeastar.f
FFILES2= table2.f makemeastar.f
FFILES3= fig25.f makemeastar.f
FOPTS  = #-w90 -w95 -w -mp -prec_div -fp_port -pc80

SRSFILES = $(FFILES)
SRSFILES2 = $(FFILES2)
SRSFILES3 = $(FFILES3)

OBJFILES = $(FFILES:.f=.o)
OBJFILES2 = $(FFILES2:.f=.o)
OBJFILES3 = $(FFILES3:.f=.o)

$(NAME):$(OBJFILES)
	$(FC) $(FOPTS) -o $(NAME) $(OBJFILES)
$(NAME2):$(OBJFILES2)
	$(FC) $(FOPTS2) -o $(NAME2) $(OBJFILES2)
$(NAME3):$(OBJFILES3)
	$(FC) $(FOPTS3) -o $(NAME3) $(OBJFILES3)

makemeastar.o : makemeastar.f params.h
	$(FC)  -c $(FOPTS) $<
driver.o : driver.f params.h
	$(FC)  -c $(FOPTS) $<
table2.o : table2.f params.h
	$(FC)  -c $(FOPTS) $<
fig25.o : fig25.f params.h
	$(FC)  -c $(FOPTS) $<
