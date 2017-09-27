CC = gcc
CFDEV = -g -O2

# CF = -ansi -std=gnu99 -pedantic -Wshadow -Wcast-qual -Wwrite-strings -g -O3
CFLAGS   = -ansi -std=gnu99
CFLAGS_1 = -ansi -pedantic -std=gnu99 -Wall -W
CFLAGS_2 = -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -g -O2

OBJECTS_contacts = get_contacts_PDB_nucl.o read_protein_PDB.o identify.o

# CFILES_contacts:
# CFILES_contacts = get_contacts_PDB_nucl.c read_protein_PDB.c identify.c
CFILES_CON := $(wildcard src_contacts/*.c)

# CFILES:
# CFILES  = main_pullends_fullGo_dcd_nucl.c read_protein.c rforce.c force_pullends_nucl.c iteration_nucl.c update_dcd_fullGo_nucl.c ras_structure_nucl.c identify.c dcdio.c
CFILES   := $(wildcard src/*.c)

# include
INC      := -Iinclude

# objects:
%.o : %.c
	$(CC) $(CFLAGS) -c $<


# CONTACTS ---------------------------------------------------Jan 29.---Mar. 21.
contacts1:
	$(CC) $(CFLAGS_1) $(CFLAGS_2) $(CFILES_CON) $(INC) -lm -o bin/run_contacts1
contacts2:
	$(CC) $(CF) $(CFILES_CON) $(INC) -lm -o bin/run_contacts2
contactADP:
	./run_contacts1 pdbADP4B9Q.pdb
contactATP:
	./run_contacts1 pdbATP4B9Q.pdb


# RUN---------------------------------------------------------Jan 29.---Mar. 21.
main:
	$(CC) $(CFLAGS_1) $(CFLAGS_2) $(CFILES) $(INC) -lm -o bin/run_protein
cluster:
	$(CC) $(CF) $(CFILES) $(INC) -lm -o bin/run_protein
m870:
	$(CC) $(CF) $(CFILES) $(INC) -lm -o bin/run_sopnucleo
m870run:
	./run_sopnucleo hsp70atp.pdb Contact_map_intra_b_70atp 1
m870dev:
	$(CC) $(CFDEV) $(CFILES) $(INC) -lm -o test/run_sopnucleo_dev
	cd test && mkdir -p Struct_data Coord
	cd test && ./run_sopnucleo_dev hsp70atpmod.pdb Contact_map_hsp70total 175

# Contact_map_intra_sc_and_b_70atp
runADP:
	./run_protein pdbADP4B9Q.pdb Contacts_intra_inter_b_ADP4B9Q_j29 1
runATP:
	./run_protein pdbATP4B9Q.pdb Contacts_intra_inter_b_ATP4B9Q_j29 1
clean:
	rm -r *.o run_protein* Struct_data/* Coord/* contact_record_*.dat out*_*.dat Info_* dihedral_potential_[0-9].dat dh_angle.dat core job-* soptime.dat *.OUT
pure:
	rm -r *.o run_protein* Struct_data/* Coord/* contact_record_*.dat out*_*.dat Info_* dihedral_potential_[0-9].dat dh_angle.dat core job* soptime.dat *.OUT
