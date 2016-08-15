#ifndef STRUCTURES_H
#define STRUCTURES_H
//Added: double occ and double beta to CAcoord;
typedef struct
{
  /*Purpose: This structure keeps CA information of a residue*/
  int    resNo;
  int    selTag; //This tag will be used to select a subset of CA coordinates. 
                 //Example: CA of a domain in protein. 
  char   residname[4];
  char   chain;
  double x;
  double y;
  double z;
  double occ;
  double beta;
}CAcoord;

typedef struct
{
  /*Purpose: This structure keeps all information of a residue*/
  int    no;
  char   chain;
  char   name[4];
  double mass;
  double charge;
  double x; //Coordinates may be coordinates of CA or Center of Mass
  double y; //Coordinates may be coordinates of CA or Center of Mass
  double z; //Coordinates may be coordinates of CA or Center of Mass
}residueInfo;

typedef struct
{
  /*Purpose: To keep all 'ATOM' records data in a pdb file according to PDB FORMAT Version 2.3*/
  int    serial;                         //Atom serial number
  char   name[5];                        //Atom name
  char   altLoc;                         //Alternate location indicator
  char   resName[4];                     //Residue name: 
  char   chainID;                        //Chain identifier
  int    resSeq;                         //Residue sequence number
  char   iCode;                          //Code for insertion of residues
  double x;                              //Orthogonal coordinates for X in Angstroms
  double y;                              //Orthogonal coordinates for Y in Angstroms
  double z;                              //Orthogonal coordinates for Z in Angstroms
  double occupancy;                      //Occupancy 
  double tempFactor;                     //Temperature factor
  char   element[2];                     //Element symbol, right-justified
  char   charge[2];                      //Charge on the atom. I dont know why it is char?????
  int    selTag;
}pdb_v23;


typedef struct
{
  /*Purpose: To keep all 'ATOM' records data in a pdb file*/
  int    resno;
  int    atnum;
  char   text[6];
  char   atom[4];
  char   attyp[5];
  char   resname[4];
  char   chain;
  char   element;
  double x;
  double y;
  double z;
  double occ;
  double temp;
}pdb2;

typedef struct
{
  /*Purpose: To keep all 'ATOM' records data in a pdb file*/
  int    resno;
  int    atnum;
  char   text[6];
  char   atom[4];
  char   attyp[5];
  char   resname[4];
  char   chain;
  char   element;
  double coordX;
  double coordY;
  double coordZ;
  double occ;
  double temp;
}pdb;


typedef struct 
{
  /*Purpose: To keep ENM residue pair information*/
  int    posI;    /* Position index of atom(residue) i.            */
  int    posJ;    /* Position index of atom(residue) j.            */
  double k;       /* Force constant between residue i and j.       */
  double dis;     /* Distance between residue i and j.             */
  double Pweight; /* A different weight can be used for different  */
		  /*  types of pairs. For example, if it contains  */
                  /*  water one can treat water pairs differently. */
}pairInfoNew;

typedef struct
{
  /*Purpose: To keep (p)air (d)istribution (f)unction of target protein*/
  double binIndex;
  double binValue; 
}pdf;

typedef struct
{
  /*Purpose: To keep energy values and then use them in Newton-Raphson method as exit condition.*/
  double energy;
}energyCutoff;

typedef struct
{
  /*Purpose: To keep secondary structure helix information.                   */
  /***********Format of Helix Information***************************************
   //Source: PDB FORMAT Version 2.3, Secondary Structure Section
   COLUMNS       DATA TYPE        FIELD        DEFINITION
   --------------------------------------------------------------------
    1 -  6       Record name      "HELIX "
    8 - 10       Integer          serNum       Serial number of the helix.
                                               This starts at 1 and increases
                                               incrementally.
   12 - 14       LString(3)       helixID      Helix identifier. In addition
                                               to a serial number, each helix is
                                               given an alphanumeric character
                                               helix identifier.
   16 - 18       Residue name     initResName  Name of the initial residue.
   20            Character        initChainID  Chain identifier for the chain
                                               containing this helix.
   22 - 25       Integer          initSeqNum   Sequence number of the initial
                                               residue.
   26            AChar            initICode    Insertion code of the initial
                                               residue.
   28 - 30       Residue name     endResName   Name of the terminal residue of
                                               the helix.
   32            Character        endChainID   Chain identifier for the chain
                                               containing this helix.
   34 - 37       Integer          endSeqNum    Sequence number of the terminal
                                               residue.
   38            AChar            endICode     Insertion code of the terminal
                                               residue.
   39 - 40       Integer          helixClass   Helix class (see below).
   41 - 70       String           comment      Comment about this helix.
   72 - 76       Integer          length       Length of this helix.
  *****************************************************************************/
  int   serNum;
  char  helixID[4];

  char  initResName[4];
  char  initChainID[2];
  int   initSeqNum;
  char  initICode[2];

  char  endResName[4];
  char  endChainID[2];
  int   endSeqNum;
  char  endICode[2];

  int   helixClass;
  /*****************************************************************************
  Helices are classified as follows:
                TYPE OF HELIX          CLASS NUMBER 
		(COLUMNS 39 - 40)
      ---------------------------------------------------
      Right-handed alpha (default)       1
      Right-handed omega                 2
      Right-handed pi                    3
      Right-handed gamma                 4
      Right-handed 310                   5
      Left-handed alpha                  6
      Left-handed omega                  7
      Left-handed gamma                  8
      27 ribbon/helix                    9
      Polyproline                       10
  *****************************************************************************/
  char  comment[31];
  int   length;
}helixInfo;

typedef struct
{
  /*Purpose: To keep secondary structure sheet information.                   */
  /***************Format of Sheet Information************************************
   //Source: PDB FORMAT Version 2.3, Secondary Structure Section
   COLUMNS     DATA TYPE        FIELD           DEFINITION
   --------------------------------------------------------------
    1 -  6     Record name      "SHEET "
    8 - 10     Integer          strand       Strand number which starts at 1 
                                             for each strand within a sheet 
                                             and increases by one.
   12 - 14     LString(3)       sheetID      Sheet identifier.
   15 - 16     Integer          numStrands   Number of strands in sheet.
   18 - 20     Residue name     initResName  Residue name of initial residue.
   22          Character        initChainID  Chain identifier of initial 
                                             residue in strand.
   23 - 26     Integer          initSeqNum   Sequence number of initial 
                                             residue in strand.
   27          AChar            initICode    Insertion code of initial residue
                                             in strand.
   29 - 31     Residue name     endResName   Residue name of terminal residue.
   33          Character        endChainID   Chain identifier of terminal
                                             residue.
   34 - 37     Integer          endSeqNum    Sequence number of terminal
                                             residue.
   38          AChar            endICode     Insertion code of terminal 
                                             residue.
   39 - 40     Integer          sense        Sense of strand with respect to
                                             previous strand in the sheet. 0
                                             if first strand, 1 if parallel,
                                             -1 if anti-parallel.
   42 - 45     Atom             curAtom      Registration. Atom name in 
                                             current strand.
   46 - 48     Residue name     curResName   Registration. Residue name in
                                             current strand.
   50          Character        curChainId   Registration. Chain identifier in
                                             current strand.
   51 - 54     Integer          curResSeq    Registration. Residue sequence
                                             number in current strand.
   55          AChar            curICode     Registration. Insertion code in
                                             current strand.
   57 - 60     Atom             prevAtom     Registration. Atom name in
                                             previous strand.
   61 - 63     Residue name     prevResName  Registration. Residue name in
                                             previous strand.
   65          Character        prevChainId  Registration. Chain identifier in
                                             previous strand.
   66 - 69     Integer          prevResSeq   Registration. Residue sequence
                                             number in previous strand.
   70          AChar            prevICode    Registration. Insertion code in
                                             previous strand.
  ******************************************************************************/
  int  strand;
  char sheetID[4];
  int  numStrands;

  char initResName[4];
  char initChainID[2];
  int  initSeqNum;
  char initICode[2];

  char endResName[4];
  char endChainID[2];
  int  endSeqNum;
  char endICode[2];

  int  sense;

  char curAtom[5];
  char curResName[4];
  char curChainID[2];
  int  curResSeq;
  char curICode[2];

  char prevAtom[5];
  char prevResName[4];
  char prevChainID[2];
  int  prevResSeq;
  char prevICode[2];
}sheetInfo;

#endif
