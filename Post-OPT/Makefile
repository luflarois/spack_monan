# Makefile para compilar post.x a partir de chem1_list.f90 e post.f90
# Compilador e flags
FC = gfortran
FFLAGS = -O2 -Wall -fbacktrace -fallow-argument-mismatch
# Para debug, descomente a linha abaixo:
# FFLAGS = -g -Wall -fbounds-check

# Nome do executável
TARGET = post.x

# Arquivos fonte
SOURCES = chem1_list.f90 post.f90
# Arquivos objeto (substitui .f90 por .o)
OBJECTS = $(SOURCES:.f90=.o)

# Regra padrão
all: $(TARGET)

# Como gerar o executável: linkando os objetos
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

# Regra implícita para .f90 -> .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Limpeza dos arquivos gerados
clean:
	rm -f $(OBJECTS) $(TARGET) *.mod

# Limpeza completa (incluindo módulos)
distclean: clean
	rm -f *.mod

# Declara que 'clean' e 'distclean' não são arquivos
.PHONY: all clean distclean
