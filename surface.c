/* #define DEBUG */
/*************************************************************************

   Program:    surface
   File:       surface.c
   
   Version:    V1.1
   Date:       19.05.94
   Function:   To create a distance map of surface features
   
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
   work! The code may not be sold commercially without prior permission 
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
   This is a very crude simple and slow version; the surface searching
   takes for ever! This could be improved dramatically by using IndexPDB
   to create an array of pointers into the linked list then sorting
   this array on each of x, y and z. Flagging the surface atoms would 
   then be extremely fast.
   
**************************************************************************

   Revision History:
   =================
   V1.0  22.11.93 Original
   V1.1  19.05.94 Added Phosphate for DNA

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines
*/
#define BOXSIZE 1.0
#define GRID    1.0
#define WATER   1.4
#define WATERSQ (WATER * WATER)

#ifdef DEBUG
#define D(BUG) fprintf(stderr,BUG)
#else
#define D(BUG)
#endif

/************************************************************************/
/* Globals
*/

/*************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL Initialise(void);
BOOL ReadCmdLine(int argc, char **argv, char *pdbfile, char *limitfile, 
                 BOOL *DoSurface);
PDB *OpenAndReadPDB(char *file, FILE *Msgfp);
PDB *FindSurfaceAtoms(PDB *pdb);
PDB *FindAtomsOfInterest(PDB *surface);
void DoDistMatrix(PDB *interest);
void Usage(void);
PDB *SelectRanges(PDB *pdb, char *limitfile);

/*************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for creating distance matrix of surface charged residues.

   18.11.93 Original   By: ACRM
   19.11.93 Modified for surface flag
*/
int main(int argc, char **argv)
{
   PDB  *pdb,
        *surface,
        *surf,
        *interest;
   char pdbfile[160],
        limitfile[160];
   BOOL DoSurface;

   if(Initialise())
   {
      if(ReadCmdLine(argc, argv, pdbfile, limitfile, &DoSurface))
      {
         if((pdb = OpenAndReadPDB(pdbfile,stderr)) != NULL)
	 {
            if(DoSurface) surface = FindSurfaceAtoms(pdb);
            else          surface = pdb;

            if(surface != NULL)
	    {
               if(limitfile[0])
                  surf = SelectRanges(surface,limitfile);
               else
                  surf = surface;

               if((interest = FindAtomsOfInterest(surf)) != NULL)
	       {
#ifdef DEBUG
                  fprintf(stderr,"\n\Interesting atom list\n");
                  WritePDB(stderr,interest);
#endif
                  DoDistMatrix(interest);
	       }
	    }
         }
      }
      else
      {
         Usage();
      }
   }

   return(0);
}

/*************************************************************************/
BOOL Initialise(void)
{
   return(TRUE);
}

/*************************************************************************/
/*>BOOL ReadCmdLine(int argc, char **argv, char *pdbfile, char *limitfile,
                    BOOL *DoSurface)
   -----------------------------------------------------------------------
   Read the command line

   18.11.93 Original   By: ACRM
   19.11.93 Added DoSurface parameter and flag
   14.06.23 Added -h for help
*/
BOOL ReadCmdLine(int argc, char **argv, char *pdbfile, char *limitfile,
                 BOOL *DoSurface)
{
   limitfile[0] = '\0';

   *DoSurface = TRUE;

   argc--; argv++;

   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
	 case 'l': case 'L':
            argc--; argv++;
            strcpy(limitfile,argv[0]);
            break;
	 case 's': case 'S':
            *DoSurface = FALSE;
            break;
         case 'h':
	 default:
            fprintf(stderr, "Returning false!\n");
            return(FALSE);
         }
      }
      else
      {
         break;
      }
      
      argc--; argv++;
   }

   if(argc != 1) return(FALSE);
   strcpy(pdbfile,argv[0]);

   return(TRUE);
}

/************************************************************************/
/*>PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
   --------------------------------------------
   Opens a PDB file specified by name and reads it with ReadPDB. Error
   messages are sent to the specified file. If NULL, no messages will
   be issued. Returns a pointer to the PDB linked list or NULL if failed.

   10.10.93 Original    By: ACRM
   11.10.93 Added legal NULL Msgfp
*/
PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
{
   FILE *fp = NULL;
   PDB  *pdb = NULL;
   int  natoms;

   if((fp=fopen(file,"r"))==NULL)
   {
      if(Msgfp!=NULL) 
         fprintf(Msgfp,"Unable to open file: %s\n",file);
   }
   else
   {
      pdb = ReadPDB(fp, &natoms);
      if(pdb == NULL || natoms == 0)
      {
         if(Msgfp!=NULL) 
            fprintf(Msgfp,"No atoms read from file: %s\n",file);
         if(pdb!=NULL) FREELIST(pdb,PDB);
         pdb = NULL;
      }
      fclose(fp);
   }
   return(pdb);
}

/*************************************************************************/
/*>PDB *FindSurfaceAtoms(PDB *pdb)
   -------------------------------
   Identifies surface atoms using a simple grid search along x, y and z
   axes. This is not an ideal analytical answer since it will not handle
   re-entrant surfaces correctly.
   Takes a PDB linked list as input and outputs a linked list of surface
   atoms.
   N.B. The occ field in the PDB linked lists will no longer be valid
   since it is used as a flag by this routine.
   18.11.93 Original   By: ACRM
*/
PDB *FindSurfaceAtoms(PDB *pdb)
{
   PDB  *surface = NULL,
        *p,
        *q;
   REAL xmin, xmax,
        ymin, ymax,
        zmin, zmax,
        x,
        y,
        z;

   fprintf(stderr,"Finding surface atoms using %f Angstrom grid\n",
           (double)GRID);

   /* Find size of coordinate box and set occ's to zero to use as flag  */
   xmin = xmax = pdb->x;
   ymin = ymax = pdb->y;
   zmin = zmax = pdb->z;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->occ = (REAL)0.0;
      if(p->x < xmin) xmin = p->x;
      if(p->x > xmax) xmax = p->x;
      if(p->y < ymin) ymin = p->y;
      if(p->y > ymax) ymax = p->y;
      if(p->z < zmin) zmin = p->z;
      if(p->z > zmax) zmax = p->z;
   }

   /* Modify these bounds by the box border size                         */
   xmin -= BOXSIZE;
   ymin -= BOXSIZE;
   zmin -= BOXSIZE;
   xmax += BOXSIZE;
   ymax += BOXSIZE;
   zmax += BOXSIZE;

   fprintf(stderr,"Dimensions of box are: %f %f %f\n\n",
           (double)(xmax-xmin),(double)(ymax-ymin),(double)(zmax-zmin));

   /* For each point on the x-y plane, search along the z axis           */
   fprintf(stderr,"Searching along z...\n");
   for(x=xmin; x<=xmax; x+=GRID)
   {
      for(y=ymin; y<=ymax; y+=GRID)
      {
         PDB grid;

         grid.x = x;
         grid.y = y;

         /* Forwards on z                                                */
         for(z=zmin; z<=zmax; z+=GRID)
         {
            BOOL GotHit=FALSE;

            grid.z = z;
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(DISTSQ(&grid, p) < WATERSQ)
 	       {
                  p->occ = (REAL)1.0;
                  GotHit = TRUE;
	       }
	    }
            if(GotHit) break;
	 }

         /* Backwards on z                                               */
         for(z=zmax; z>=zmin; z-=GRID)
         {
            BOOL GotHit=FALSE;

            grid.z = z;
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(DISTSQ(&grid, p) < WATERSQ)
 	       {
                  p->occ = (REAL)1.0;
                  GotHit = TRUE;
	       }
	    }
            if(GotHit) break;
	 }
      }
   }


   /* For each point on the x-z plane, search along the y axis           */
   fprintf(stderr,"Searching along y...\n");
   for(x=xmin; x<=xmax; x+=GRID)
   {
      for(z=zmin; z<=zmax; z+=GRID)
      {
         PDB grid;

         grid.x = x;
         grid.z = z;

         /* Forwards on y                                                */
         for(y=ymin; y<=ymax; y+=GRID)
         {
            BOOL GotHit=FALSE;

            grid.y = y;
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(DISTSQ(&grid, p) < WATERSQ)
 	       {
                  p->occ = (REAL)1.0;
                  GotHit = TRUE;
	       }
	    }
            if(GotHit) break;
	 }

         /* Backwards on y                                               */
         for(y=ymax; y>=zmin; y-=GRID)
         {
            BOOL GotHit=FALSE;

            grid.y = y;
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(DISTSQ(&grid, p) < WATERSQ)
 	       {
                  p->occ = (REAL)1.0;
                  GotHit = TRUE;
	       }
	    }
            if(GotHit) break;
	 }
      }
   }


   /* For each point on the y-z plane, search along the x axis           */
   fprintf(stderr,"Searching along x...\n");
   for(y=ymin; y<=ymax; y+=GRID)
   {
      for(z=zmin; z<=zmax; z+=GRID)
      {
         PDB grid;

         grid.z = z;
         grid.y = y;

         /* Forwards on x                                                */
         for(x=xmin; x<=xmax; x+=GRID)
         {
            BOOL GotHit=FALSE;

            grid.x = x;
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(DISTSQ(&grid, p) < WATERSQ)
 	       {
                  p->occ = (REAL)1.0;
                  GotHit = TRUE;
	       }
	    }
            if(GotHit) break;
	 }

         /* Backwards on x                                               */
         for(x=xmax; x>=xmin; x-=GRID)
         {
            BOOL GotHit=FALSE;

            grid.x = x;
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(DISTSQ(&grid, p) < WATERSQ)
 	       {
                  p->occ = (REAL)1.0;
                  GotHit = TRUE;
	       }
	    }
            if(GotHit) break;
	 }
      }
   }

   /* Now copy the flagged atoms into an output linked list              */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->occ != (REAL)0.0)
      {
         if(surface == NULL)
         {
            INIT(surface,PDB);
            q=surface;
         }
         else
         {
            ALLOCNEXT(q,PDB);
         }

         if(q==NULL)
	 {
            if(surface!=NULL) FREELIST(surface,PDB);
            fprintf(stderr,"No memory for surface list\n");
            return(NULL);
         }

         CopyPDB(q,p);
      }
   }

   /* And return the linked list of flagged atoms                        */
   return(surface);
}

/*************************************************************************/
/*>PDB *FindAtomsOfInterest(PDB *surface)
   --------------------------------------
   Searches a PDB linked list for all charged atoms and returns a PDB 
   linked list containing only those atoms.
   This could be improved to read the atoms of interest from a file
   for more flexibility. At the moment, we just stick with charged
   atoms.

   18.11.93 Original   By: ACRM
   22.11.93 Added aromatics
   19.05.94 Added phosphate for DNA
*/
PDB *FindAtomsOfInterest(PDB *surface)
{
   PDB *interest = NULL,
       *CurrInt  = NULL,
       *p,
       *q;

   fprintf(stderr,"Finding atoms of interest on the surface...\n");

   /* We use occ as a flag to indicate interesting atoms                 */
   D("Clearing occ flag\n");
   for(p=surface; p!=NULL; NEXT(p))
      p->occ = 0.0;

   D("Finding charged and aromatic residues\n");
   for(p=surface; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->resnam,"GLU",3) && !strncmp(p->atnam,"OE",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"ASP",3) && !strncmp(p->atnam,"OD",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"ARG",3) && !strncmp(p->atnam,"NE",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"ARG",3) && !strncmp(p->atnam,"CZ",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"ARG",3) && !strncmp(p->atnam,"NH",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"LYS",3) && !strncmp(p->atnam,"NZ",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"PHE",3) && !strncmp(p->atnam,"CG",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"PHE",3) && !strncmp(p->atnam,"CD",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"PHE",3) && !strncmp(p->atnam,"CE",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"PHE",3) && !strncmp(p->atnam,"CZ",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TYR",3) && !strncmp(p->atnam,"CG",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TYR",3) && !strncmp(p->atnam,"CD",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TYR",3) && !strncmp(p->atnam,"CE",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TYR",3) && !strncmp(p->atnam,"CZ",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TRP",3) && !strncmp(p->atnam,"CD",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TRP",3) && !strncmp(p->atnam,"NE",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TRP",3) && !strncmp(p->atnam,"CE",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TRP",3) && !strncmp(p->atnam,"CZ",2))
         p->occ = 1.0;
      if(!strncmp(p->resnam,"TRP",3) && !strncmp(p->atnam,"CH",2))
         p->occ = 1.0;
      if(!strncmp(p->atnam,"P ",2))
         p->occ = 1.0;
   }

   /* Now we copy this list, but include only one entry for each residue */
   for(p=surface; p!=NULL; NEXT(p))
   {
      if(p->occ != 0.0)
      {
         BOOL ResFound = FALSE;

         /* See if this residue has been flagged already                 */
         if(interest != NULL)
         {
            D("Searching for this residue found already\n");
            for(q=interest; q!=NULL; NEXT(q))
            {
               if(q->resnum    == p->resnum    && 
                  q->insert[0] == p->insert[0] &&
                  q->chain[0]  == p->chain[0])
               {
                  /* Yes, it's been found before so shift CofG           */
                  D("   YES: updating CofG\n");

                  q->x *= q->occ;
                  q->x += p->x;
                  (q->occ) += (REAL)1.0;
                  q->x /= q->occ;

                  q->y *= q->occ;
                  q->y += p->y;
                  (q->occ) += (REAL)1.0;
                  q->y /= q->occ;

                  q->z *= q->occ;
                  q->z += p->z;
                  (q->occ) += (REAL)1.0;
                  q->z /= q->occ;

                  ResFound = TRUE;
                  break;
               }
            }
         }

         /* If not found, then add to the list                          */
         if(!ResFound)
	 {
            D("   Not Found: adding to list\n");

            if(interest == NULL)
            {
               INIT(interest,PDB);
               CurrInt = interest;
            }
            else
            {
               ALLOCNEXT(CurrInt,PDB);
            }

            if(CurrInt==NULL)
	    {
               if(interest!=NULL) FREELIST(interest,PDB);
               fprintf(stderr,"No memory for charged atom list\n");
               return(NULL);
            }

            CopyPDB(CurrInt,p);
            /* Use occ to store a count of number of times this residue
               has been found.
            */
            CurrInt->occ = (REAL)1.0;
	 }
      }
   }
   return(interest);
}

/*************************************************************************/
/*>void DoDistMatrix(PDB *interest)
   --------------------------------
   Calculate a distance matrix for all the atoms in the input PDB linked
   list. Print this to stdout.

   18.11.93 Original   By: ACRM
*/
void DoDistMatrix(PDB *interest)
{
   PDB *p, *q;

   D("Entered DoDistMatrix\n");

   for(p=interest; p!=NULL; NEXT(p))
   {
      for(q=p->next; q!=NULL; NEXT(q))
      {
         printf("%s %c%4d%c %s %c%4d%c %8.3f\n",
                p->resnam,p->chain[0],p->resnum,p->insert[0],
                q->resnam,q->chain[0],q->resnum,q->insert[0],
                DIST(p,q));
      }
   }
}

/*************************************************************************/
/*>void Usage(void)
   ----------------
   Display a usage message

   18.11.93 Original   By: ACRM
   19.11.93 Added -s flag
*/
void Usage(void)
{
   fprintf(stderr,"Surface V1.1 19.05.93\n");

   fprintf(stderr,"\nUsage: surface [-l <limits>] [-s] <file.pdb>\n");
   fprintf(stderr,"       -l specify limits file\n");
   fprintf(stderr,"       -s assume all residues are surface\n");
   fprintf(stderr,"\nSearch for surface charged and aromatic residues and \
create a distance matrix\n");
   fprintf(stderr,"The limits file specifies ranges of amino acids to \
be included\n");
   fprintf(stderr,"Output is to stdout\n");
}

/*************************************************************************/
/*>PDB *SelectRanges(PDB *pdb, char *limitfile)
   --------------------------------------------
   Read the file specified by limitfile containing amino acid residue
   ranges specified as:
      [c]nnn[i]  [c]nnn[i]
   where [c] is an optional chain specification, nnn is a residue number
   and [i] is an optional insert specification.
   Then create a new PDB linked list from the input list, containing only
   those residues which fall in the specified ranges.

   18.11.93 Original   By: ACRM
*/
PDB *SelectRanges(PDB *pdb, char *limitfile)
{
   FILE *fp;
   PDB  *start,
        *stop,
        *p,
        *q,
        *outpdb = NULL;
   char buffer[160],
        spec1[16],
        spec2[16];
   int  nrange=0,
        j;

   /* First read the range limits file                                   */
   if((fp=fopen(limitfile,"r"))==NULL)
   {
      fprintf(stderr,"Unable to read limits file (ignored).\n");
      return(pdb);
   }

   /* Scan through to count ranges                                       */
   while(fgets(buffer,160,fp)) nrange++;
   rewind(fp);

   /* Allocate memory for the range arrays                               */
   if((start = (PDB *)malloc(nrange * sizeof(PDB))) == NULL)
   {
      fprintf(stderr,"No memory for limit range storage. Limits ignored\n");
      return(pdb);
   }   
   if((stop = (PDB *)malloc(nrange * sizeof(PDB))) == NULL)
   {
      free(stop);
      fprintf(stderr,"No memory for limit range storage. Limits ignored\n");
      return(pdb);
   }

   /* Read the file again, storing ranges                                */
   nrange = 0;
   while(fgets(buffer,160,fp))
   {
      if(sscanf(buffer,"%s %s",spec1, spec2) == 2)
      {
         char chain[8],
              insert[8];
         int  resnum;

         ParseResSpec(spec1,chain,&resnum,insert);
         start[nrange].chain[0]  = chain[0];
         start[nrange].chain[1]  = '\0';
         start[nrange].insert[0] = insert[0];
         start[nrange].insert[1] = '\0';
         start[nrange].resnum    = resnum;

         ParseResSpec(spec2,chain,&resnum,insert);
         stop[nrange].chain[0]   = chain[0];
         stop[nrange].chain[1]   = '\0';
         stop[nrange].insert[0]  = insert[0];
         stop[nrange].insert[1]  = '\0';
         stop[nrange].resnum     = resnum;

         nrange++;
      }
   }

   fclose(fp);

   /* Now we've read the ranges, create a new PDB linked list containing
      only residues which fall inside these ranges
   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* See if this record is in any of the ranges                      */
      for(j=0; j<nrange; j++)
      {
         if(p->chain[0] == start[j].chain[0] &&
            p->resnum   >= start[j].resnum   &&
            p->resnum   <= stop[j].resnum)
	 {
            /* Numbers and chain are OK, check for insert                */
            if(p->resnum == start[j].resnum  &&
               start[j].insert[0] != ' '     &&
               p->insert[0] < start[j].insert[0])
               continue;
            if(p->resnum == stop[j].resnum   &&
               stop[j].insert[0] != ' '      &&
               p->insert[0] > stop[j].insert[0])
               continue;

            /* This one is within ranges, so create space and copy       */
            if(outpdb==NULL)
	    {
               INIT(outpdb,PDB);
               q=outpdb;
	    }
            else
	    {
               ALLOCNEXT(q,PDB);
	    }

            if(q==NULL)
	      {
               if(outpdb!=NULL) FREELIST(outpdb,PDB);
               fprintf(stderr,"No memory for residues in specified zones. \
Using all residues.\n");
               return(pdb);
	      }

            CopyPDB(q,p);

            /* Break out of search through ranges                        */
            break;
         }
      }
   }

   /* Free memory and return                                             */
   if(start != NULL) free(start);
   if(stop  != NULL) free(stop);

   return(outpdb);
}
