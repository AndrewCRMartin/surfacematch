/*************************************************************************

   Program:    match
   File:       match.c
   
   Version:    V1.0
   Date:       21.11.93
   Function:   Match 2 distance matrices as created by surface
   
   Copyright:  (c) SciTech Software 1993
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is released under the GPL3 licence with the following
   provisos:

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without print permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and that 
   the source code is provided with the program.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  21.11.93 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/fsscanf.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines
*/
#define MAXPROP         8
#define PROP_POSITIVE   0
#define PROP_NEGATIVE   1
#define PROP_AROM       2

#define MAXDIST        32

#define DEFBIN        1.0     /* Default distance bin size              */
#define DEFACC       50.0     /* Default string match accuracy          */

/************************************************************************/
/* Structure and type definitions
*/
typedef struct
{
   int  resnum[2];
   char resnam[2][8],
        chain[2][8],
        insert[2][8];
   int  dist;
   BOOL dead;
}  DATA;

typedef struct
{
   int  resnum;
   char resnam[8],
        chain[8],
        insert[8];
   char dist[MAXDIST],
        property[MAXPROP];
}  ATOM;

/*************************************************************************/
/* Globals
*/
REAL gBin      = DEFBIN,    /* Bin size for distance matrix              */
     gAccuracy = DEFACC;    /* Percentage accuracy for string comparison */

/*************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ReadCmdLine(int argc, char **argv, char *PatFile, char *StrucFile);
void Usage(void);
void MatchFiles(FILE *fp_pat, FILE *fp_struc);
DATA *ReadMatrix(FILE *fp, int *outnatom);
ATOM *CreateAtomArray(DATA *data, int ndata, int *outnatom, BOOL SwapProp);
int ConvDist(REAL dist);
int GotAtom(ATOM *outdata, int natom, char *chain, int resnum, 
            char *insert);
void FillAtom(ATOM *outdata, int natom, int pos, char *chain, int resnum, 
              char *insert, char *resnam, int DistRange, BOOL SwapProp);
BOOL DoLesk(int npat, DATA *pat, int nstruc, DATA *struc);
void KillAtom(char *chain, int resnum, char *insert, DATA *data, 
              int ndata);
void TrimBitStrings(int npat, ATOM *pat, int nstruc, ATOM *struc);
void PrintResults(int NPatAtom,   ATOM *PatAtom, 
                  int NStrucAtom, ATOM *StrucAtom);
void PrintBestMatch(ATOM *PatAtom,   int PatIndex, 
                    ATOM *StrucAtom, int NStrucAtom);
BOOL Compare(char *str1, char *str2, int len);
REAL CalcScore(char *str1, char *str2, int len);

/*************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for matching output files from surface.

   18.11.93 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char PatFile[160],
        StrucFile[160];
   FILE *fp_pat,
        *fp_struc;

   if(ReadCmdLine(argc, argv, PatFile, StrucFile))
   {
      if((fp_pat = fopen(PatFile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open pattern file: %s\n",PatFile);
         exit(1);
      }
      if((fp_struc = fopen(StrucFile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open pattern file: %s\n",StrucFile);
         exit(1);
      }

      MatchFiles(fp_pat, fp_struc);
   }
   else
   {
      Usage();
   }
   return(0);
}

/*************************************************************************/
/*>BOOL ReadCmdLine(int argc, char **argv, char *PatFile, char *StrucFile)
   -----------------------------------------------------------------------
   Read the command line

   18.11.93 Original   By: ACRM
   22.11.93 Added flags
*/
BOOL ReadCmdLine(int argc, char **argv, char *PatFile, char *StrucFile)
{
   argc--; argv++;

   while(argc>2)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd': case 'D':
            argc--; argv++;
            sscanf(argv[0],"%lf",&gBin);
            if(gBin == 0.0) gBin = 1.0;
            break;
         case 'a': case 'A':
            argc--; argv++;
            sscanf(argv[0],"%lf",&gAccuracy);
            if(gAccuracy == 0.0) gAccuracy = 100.0;
            break;
	 default:
            break;
         }
      }
      argc--; argv++;
   }

   if(argc != 2) return(FALSE);

   strcpy(PatFile,argv[0]);
   strcpy(StrucFile,argv[1]);

   return(TRUE);
}

/*************************************************************************/
/*>void Usage(void)
   ----------------
   Display a usage message

   18.11.93 Original   By: ACRM
   22.11.93 Added flag decriptions
*/
void Usage(void)
{
   fprintf(stderr,"Usage: match [-d <bin>] [-a <accuracy>] <pattern file> \
<structure file>\n");
   fprintf(stderr,"       -d specifies distance bin size (default: %.1lf)\n",
           (double)DEFBIN);
   fprintf(stderr,"       -a specifies percent dist string match accuracy \
(default: %.1lf)\n",(double)DEFACC);
   fprintf(stderr,"\nMatch V1.0 21.11.93\n");
   fprintf(stderr,"Find potential matches for a pattern in a structure \
using Lesks's method\n");
   fprintf(stderr,"The input files are generated by Surface V1.0\n");
   fprintf(stderr,"Currently handles properties for charged and aromatic \
residues\n");
}

/*************************************************************************/
/*>void MatchFiles(FILE *fp_pat, FILE *fp_struc)
   ---------------------------------------------
   18.11.93 Original   By: ACRM
*/
void MatchFiles(FILE *fp_pat, FILE *fp_struc)
{
   DATA *pat,
        *struc;
   int  npat,
        nstruc;

   pat   = ReadMatrix(fp_pat,   &npat);
   struc = ReadMatrix(fp_struc, &nstruc);

   printf("%d records read from pattern; %d records read from structure\n",
          npat, nstruc);

   DoLesk(npat, pat, nstruc, struc);
}

/*************************************************************************/
/*>ATOM *ReadMatrix(FILE *fp, int *natom)
   --------------------------------------
   Read the output from surface. Create an array of type ATOM which
   contains flags for the presence of connections of a given distance and
   flags the atom type.

   18.11.93 Original   By: ACRM
   22.11.93 Corrected return values
*/
DATA *ReadMatrix(FILE *fp, int *outnatom)
{
   int  maxrec   = 0,
        i        = 0;
   DATA *outdata = NULL;
   char buffer[160];

   *outnatom = 0;

   /* Scan the file to see how long it is                                */
   while(fgets(buffer,160,fp)) maxrec++;
   rewind(fp);

   /* Allocate this much space                                           */
   if((outdata = malloc(maxrec * sizeof(DATA)))==NULL)
      return(NULL);

   /* Now read the file again, filling in data                           */
   i=0;
   while(fgets(buffer,160,fp))
   {
      char resnam1[8], resnam2[8],
           chain1[8],  chain2[8],
           insert1[8], insert2[8];
      int  resnum1,    resnum2,
           DistRange;
      REAL dist;

      fsscanf(buffer,"%4s%1x%c%4d%c%1x%4s%1x%c%4d%c%1x%8lf",
              resnam1, chain1, &resnum1, insert1,
              resnam2, chain2, &resnum2, insert2,
              &dist);

      DistRange = ConvDist(dist);

      outdata[i].resnum[0] = resnum1;
      strcpy(outdata[i].resnam[0], resnam1);
      strcpy(outdata[i].chain[0],  chain1);
      strcpy(outdata[i].insert[0], insert1);

      outdata[i].resnum[1] = resnum2;
      strcpy(outdata[i].resnam[1], resnam2);
      strcpy(outdata[i].chain[1],  chain2);
      strcpy(outdata[i].insert[1], insert2);

      outdata[i].dist = DistRange;
      outdata[i].dead = FALSE;

      i++;
   }

   *outnatom = i;
   return(outdata);
}

/*************************************************************************/
/*>ATOM *CreateAtomArray(DATA *data, int ndata, int *outnatom, 
                         BOOL SwapProp)
   -----------------------------------------------------------
   19.11.93 Original   By: ACRM
*/
ATOM *CreateAtomArray(DATA *data, int ndata, int *outnatom, BOOL SwapProp)
{
   int  i,
        natom = 0,
        pos;
   ATOM *outatom;

   /* Allocate memory for the output atom array                          */
   if((outatom = (ATOM *)malloc(ndata * sizeof(ATOM))) == NULL)
      return(NULL);

   for(i=0; i<ndata; i++)
   {
      /* Skip this record if the dead flag is set                        */
      if(data[i].dead) continue;

      /* See if we've got a record for the first residue                 */
      pos = GotAtom(outatom, natom, data[i].chain[0], data[i].resnum[0],
                    data[i].insert[0]);
      if(pos == (-1)) pos = natom++;

      /* Fill this into the data array                                   */
      FillAtom(outatom, natom, pos, data[i].chain[0], data[i].resnum[0],
               data[i].insert[0], data[i].resnam[0], data[i].dist, SwapProp);

      /* See if we've got a record for the second residue                */
      pos = GotAtom(outatom, natom, data[i].chain[1], data[i].resnum[1],
                    data[i].insert[1]);
      if(pos == (-1)) pos = natom++;

      /* Fill this into the data array                                   */
      FillAtom(outatom, natom, pos, data[i].chain[1], data[i].resnum[1],
               data[i].insert[1], data[i].resnam[1], data[i].dist, SwapProp);
   }

   *outnatom = natom;

   return(outatom);
}

/*************************************************************************/
/*>int ConvDist(REAL dist)
   -----------------------
   Converts a distance (REAL) to an integer bin number < MAXDIST. The bin
   size os read from the global variable gBin

   19.11.93 Original   By: ACRM
   22.11.93 Changed to read bin size from global gBin
*/
int ConvDist(REAL dist)
{
   int idist;

   idist = (int)(dist/gBin);
   if(idist >= MAXDIST-1) idist = MAXDIST-2;

   return(idist);
}

/*************************************************************************/
/*>int GotAtom(ATOM *outdata, int natom, char *chain, int resnum, 
               char *insert)
   --------------------------------------------------------------
   Searches the current outdata array to see if we already have this
   residue. If so, returns the index into the array. If not, returns -1

   19.11.93 Original   By: ACRM
*/
int GotAtom(ATOM *outdata, int natom, char *chain, int resnum, 
            char *insert)
{
   int i;

   for(i=0; i<natom; i++)
   {
      if(outdata[i].resnum    == resnum    &&
         outdata[i].chain[0]  == chain[0]  &&
         outdata[i].insert[0] == insert[0])
      {
         return(i);
      }
   }
   return(-1);
}

/*************************************************************************/
/*>void FillAtom(ATOM *outdata, int natom, int pos, char *chain, 
                 int resnum, char *insert, char *resnam, int DistRange,
                 BOOL SwapProp)
   --------------------------------------------------------------------
   Fill in an item in the data array. If (pos == natom-1) then it's a
   new residue so we must fill in all data; otherwise just set the
   appropriate flags.

   19.11.93 Original   By: ACRM
   22.11.93 Added aromatic support
   19.05.94 Added DNA support
*/
void FillAtom(ATOM *outdata, int natom, int pos, char *chain, int resnum, 
              char *insert, char *resnam, int DistRange, BOOL SwapProp)
{
   int i;

   if(pos == natom-1)
   {
      /* A new residue; fill in all data                                 */
      outdata[pos].resnum = resnum;
      strcpy(outdata[pos].chain,  chain);
      strcpy(outdata[pos].insert, insert);
      strcpy(outdata[pos].resnam, resnam);

      /* Clear all distance and property flags                           */
      for(i=0; i<MAXPROP; i++) outdata[pos].property[i] = '0';
      outdata[pos].property[MAXPROP-1] = '\0';
      for(i=0; i<MAXDIST; i++) outdata[pos].dist[i]     = '0';
      outdata[pos].dist[MAXDIST-1] = '\0';

      /* Set the property flags for this residue                         */
      if(!strncmp(resnam,"LYS",3) || !strncmp(resnam,"ARG",3))
         outdata[pos].property[SwapProp?PROP_NEGATIVE:PROP_POSITIVE] = '1';
      if(!strncmp(resnam,"GLU",3) || 
         !strncmp(resnam,"ASP",3) ||
         !strncmp(resnam,"  A",3) ||
         !strncmp(resnam,"  T",3) ||
         !strncmp(resnam,"  C",3) ||
         !strncmp(resnam,"  G",3))
         outdata[pos].property[SwapProp?PROP_POSITIVE:PROP_NEGATIVE] = '1';
      if(!strncmp(resnam,"TYR",3) || 
         !strncmp(resnam,"PHE",3) ||
         !strncmp(resnam,"TRP",3))
         outdata[pos].property[PROP_AROM] = '1';
   }

   /* Now set the distance flag                                          */
   outdata[pos].dist[DistRange] = '1';
}

/*************************************************************************/
/*>BOOL DoLesk(int npat, DATA *pat, int nstruc, DATA *struc)
   ---------------------------------------------------------
   Does the actual Lesk pattern matching algorithm (with some 
   modifications).

   19.11.93 Original   By: ACRM
   21.11.93 Added property comparison and printing of results :-)
*/
BOOL DoLesk(int npat, DATA *pat, int nstruc, DATA *struc)
{
   ATOM *PatAtom,
        *StrucAtom;
   int  NPatAtom,
        NStrucAtom,
        i, j, k,
        PrevPatAtoms   = 0,
        PrevStrucAtoms = 0;
   
   for(i=0; i<100; i++)
   {
      /* Create the bit strings for the atoms from the data arrays       */
      PatAtom   = CreateAtomArray(pat,   npat,   &NPatAtom, TRUE);
      StrucAtom = CreateAtomArray(struc, nstruc, &NStrucAtom, FALSE);

      /* Print information on remaining atoms                            */
      printf("Iteration %d: %d pattern atoms and %d structure atoms \
remain\n",i,NPatAtom,NStrucAtom);

      /* Remove any distance flags from the bit strings which are 
         never seen in the other structure
      */
      TrimBitStrings(NPatAtom, PatAtom, NStrucAtom, StrucAtom);

      /* Exit if we've converged                                         */
      if(NPatAtom == PrevPatAtoms && NStrucAtom == PrevStrucAtoms)
         break;
      PrevPatAtoms   = NPatAtom;
      PrevStrucAtoms = NStrucAtom;
         
#ifdef DEBUG
      printf("\nPattern atoms are:\n");
      for(j=0;j<NPatAtom;j++)
      {
         printf("Property: %s; Distance: %s\n",
                PatAtom[j].property, PatAtom[j].dist);
      }

      printf("\nStructure atoms are:\n");
      for(j=0;j<NStrucAtom;j++)
      {
         printf("Property: %s; Distance: %s\n",
                StrucAtom[j].property, StrucAtom[j].dist);
      }
#endif      

      /* For each atom in the structure look to see if the bit string is
         not found in the pattern. If not found, kill the atom
      */
      for(j=0; j<NStrucAtom; j++)
      {
         BOOL Found = FALSE;

         for(k=0; k<NPatAtom; k++)
         {
	    if(!strncmp(StrucAtom[j].property,PatAtom[k].property,MAXPROP))
	    {
	       if(!Compare(StrucAtom[j].dist,PatAtom[k].dist,MAXDIST))
	       {
		  Found = TRUE;
		  break;
	       }
	    }
	 }

         if(!Found)
	 {
            KillAtom(StrucAtom[j].chain, StrucAtom[j].resnum, 
                     StrucAtom[j].insert, struc, nstruc);
	    printf("Structure atom: %s %c%d%c killed\n",
		   StrucAtom[j].resnam, StrucAtom[j].chain[0],
		   StrucAtom[j].resnum, StrucAtom[j].insert[0]);
         }
      }
   }

   PrintResults(NPatAtom, PatAtom, NStrucAtom, StrucAtom);
}

/*************************************************************************/
/*>void KillAtom(char *chain, int resnum, char *insert, 
                 DATA *struc, int nstruc)
   ----------------------------------------------------
   Set the dead flag in the data array for an atom;

   19.11.93 Original   By: ACRM
*/
void KillAtom(char *chain, int resnum, char *insert, DATA *data, 
              int ndata)
{
   int i;

   for(i=0; i<ndata; i++)
   {
      if((resnum    == data[i].resnum[0] &&
          chain[0]  == data[i].chain[0][0]  &&
          insert[0] == data[i].insert[0][0])   ||
         (resnum    == data[i].resnum[1] &&
          chain[0]  == data[i].chain[1][0]  &&
          insert[0] == data[i].insert[1][0]))
      {
         data[i].dead = TRUE;
      }
   }
}

/*************************************************************************/
/*>void TrimBitStrings(int npat, ATOM *pat, int nstruc, ATOM *struc)
   -----------------------------------------------------------------
   Search the bit strings of the pattern and remove any distance flags
   which never occur in the structure and vice versa

   19.11.93 Original   By: ACRM
*/
void TrimBitStrings(int npat, ATOM *pat, int nstruc, ATOM *struc)
{
   int  i, j, k;
   BOOL PatHit,
        StrucHit;

   /* For each distance in the distance bit string                       */
   for(i=0; i<MAXDIST; i++)
   {
      PatHit = StrucHit = FALSE;

      /* Search the pattern array for this distance being flagged        */
      for(j=0; j<npat; j++)
      {
         if(pat[j].dist[i] == '1')
	 {
            PatHit = TRUE;
            break;
         }
      }

      /* Search the structure array for this distance being flagged      */
      for(j=0; j<nstruc; j++)
      {
         if(struc[j].dist[i] == '1')
	 {
            StrucHit = TRUE;
            break;
         }
      }

      /* If there was a hit in pattern, but not structure, kill refs in
         pattern
      */
      if(PatHit && !StrucHit)
      {
         for(j=0; j<npat; j++)
	 {
            pat[j].dist[i] = '0';
         }
      }

      /* If there was a hit in structure, but not pattern, kill refs in
         structure
      */
      if(StrucHit && !PatHit)
      {
         for(j=0; j<nstruc; j++)
	 {
            struc[j].dist[i] = '0';
         }
      }
   }
}

/*************************************************************************/
/*>void PrintResults(int NPatAtom,   ATOM *PatAtom, 
                     int NStrucAtom, ATOM *StrucAtom)
   --------------------------------------------------
   Run through the pattern atoms and, for each, print the best match 
   from the structure atoms

   22.11.93 Original   By: ACRM
*/
void PrintResults(int NPatAtom,   ATOM *PatAtom, 
                  int NStrucAtom, ATOM *StrucAtom)
{
   int i,
       j;

   for(i=0; i<NPatAtom; i++)
      PrintBestMatch(PatAtom, i, StrucAtom, NStrucAtom);
}

/*************************************************************************/
/*>void PrintBestMatch(ATOM *PatAtom,   int PatIndex, 
                       ATOM *StrucAtom, int NStrucAtom)
   ----------------------------------------------------
   Print the best macth from the structure for this pattern atom

   22.11.93 Original   By: ACRM
*/
void PrintBestMatch(ATOM *PatAtom,   int PatIndex, 
                    ATOM *StrucAtom, int NStrucAtom)
{
   int  j,
        best = -1;
   REAL score, 
        BestScore = 0.0;

   for(j=0; j<NStrucAtom; j++)
   {
      if(!strncmp(PatAtom[PatIndex].property,StrucAtom[j].property,MAXPROP) &&
         !Compare(PatAtom[PatIndex].dist,    StrucAtom[j].dist,    MAXDIST))
      {
         score = CalcScore(PatAtom[PatIndex].dist, StrucAtom[j].dist, MAXDIST);
         if(score > BestScore)
         {
            BestScore = score;
            best = j;
         }
      }
   }

   /* If we got a best score, print it out                              */
   if(best != (-1))
   {
      printf("Pattern: %s %c%d%c matches Structure: %s %c%d%c\n",
             PatAtom[PatIndex].resnam, PatAtom[PatIndex].chain[0], 
             PatAtom[PatIndex].resnum, PatAtom[PatIndex].insert[0],
             StrucAtom[best].resnam, StrucAtom[best].chain[0], 
             StrucAtom[best].resnum, StrucAtom[best].insert[0]);
   }
}

/*************************************************************************/
/*>BOOL Compare(char *str1, chsr *str2, int len)
   ---------------------------------------------
   Compares two bitstrings (1/0 character strings) using a 75% identity
   requirement

   22.11.93 Original   By: ACRM
*/
BOOL Compare(char *str1, char *str2, int len)
{
   int  i,
        CountStr1 = 0,
        CountStr2 = 0,
        NMatch    = 0;
   REAL Accuracy;

   for(i=0; i<len; i++)
   {
      if(str1[i] == '1') CountStr1++;
      if(str2[i] == '1') CountStr2++;

      if((str1[i] == '1') && (str2[i] == '1')) NMatch++;
   }

   if((Accuracy = ((REAL)100.0 * (REAL)NMatch / 
                   (REAL)MAX(CountStr1, CountStr2))) >= 
      gAccuracy)
   {
      return(FALSE);
   }

   return(TRUE);
}

/*************************************************************************/
/*>REAL CalcScore(char *str1, char *str2, int len)
   -----------------------------------------------
   Calculate the score for this match. Much the same as compare, but 
   returns the percentage score rather than a BOOL

   22.11.93 Original   By: ACRM
*/
REAL CalcScore(char *str1, char *str2, int len)
{
   int i,
       CountStr1 = 0,
       CountStr2 = 0,
       NMatch    = 0;

   for(i=0; i<len; i++)
   {
      if(str1[i] == '1') CountStr1++;
      if(str2[i] == '1') CountStr2++;

      if((str1[i] == '1') && (str2[i] == '1')) NMatch++;
   }

   return((REAL)100.0 * (REAL)NMatch / (REAL)MAX(CountStr1, CountStr2));
}
