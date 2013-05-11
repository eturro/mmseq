/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "fasta.h"


/* Function: OpenFASTA(), ReadFASTA(), CloseFASTA().
 * Date:     SRE, Sun Sep  8 06:39:26 2002 [AA2721, transatlantic]
 *
 * Purpose:  A very rudimentary FASTA file reading API. Designed
 *           for simplicity and clarity, not for robustness.
 *           
 *           The API is:
 *           
 *           ffp = OpenFASTA(seqfile);
 *           while (ReadFASTA(ffp, &seq, &name, &seqlen)
 *           {
 *             do stuff with sequence;
 *             free(name);
 *             free(seq);
 *           }
 *           CloseFASTA(ffp);
 *           
 * Args:     
 *           seqfile   - name of a FASTA file to open.
 *           seq       - RETURN: one sequence
 *           name      - RETURN: name of the sequence
 *           seqlen    - RETURN: length of the sequence in residues
 *           ffp       - ptr to a FASTAFILE object.
 *           
 * Commentary: 
 *           The basic problem with reading FASTA files is that there is
 *           no end-of-record indicator. When you're reading sequence n, 
 *           you don't know you're done until you've read the header line
 *           for sequence n+1, which you won't parse 'til later (when 
 *           you're reading in the sequence n+1). One common trick for
 *           this is to implement a one-line "lookahead" buffer that you
 *           can peek at, before parsing later. 
 *           
 *           This buffer is kept in a small structure (a FASTAFILE), rather
 *           than in a static char[] in the function. This allows
 *           us to have multiple FASTA files open at once. The static approach
 *           would only allow us to have one file open at a time. ANSI C
 *           predates the widespread use of parallel programming. It was
 *           not overly concerned about the drawbacks of statics. Today,
 *           though, you should keep in mind that you may someday want to
 *           turn your program into a multithreaded, parallel program, and
 *           all functions in parallelized code must be "reentrant": able to
 *           be called a second time - with different arguments,
 *           and while the code in the first function call is still executing! -
 *           without overwriting or corrupting any static storage in the
 *           function. Statics have fewer uses now (for example, to
 *           test that some initialization code for a function is run once 
 *           and only once.)
 * 
 * Limitations:          
 *           There is no error handling, for clarity's sake. Also,
 *           the parser is brittle. Improper FASTA files (for instance,
 *           blank lines between records) will cause unexpected
 *           behavior. Real file parsers are more complex.
 *           In real life, they have to deal with absolutely anything the user might
 *           pass as a "FASTA file"; and either parse it correctly,
 *           or detect that it's an invalid format and fail cleanly.
 *           
 *           Lines are read in from the file using ANSI C's fgets(). fgets()
 *           requires a maximum buffer length (here, FASTA_MAXLINE, which is
 *           defined as 512 in bio5495.h). Some FASTA files have very long
 *           description lines, however; notably the NCBI NR database. Static
 *           limitations on things like line or sequence lengths should be
 *           avoided. An example of a replacement for fgets() that dynamically
 *           allocates its buffer size and allows any line length is
 *           SQUID's sre_fgets().
 *           
 *           We use ANSI C's strtok() to parse the sequence name out of the line.
 *           strtok() is deprecated in modern programs because it is not threadsafe. 
 *           (See comments above.) An example of a threadsafe version is
 *           SQUID's sre_strtok().
 *           
 * Returns:  
 *           OpenFASTA() returns a FASTAFILE pointer, or NULL on failure (for
 *           instance, if the file doesn't exist, or isn't readable).
 *           
 *           ReadFASTA() returns 1 on success, or a 0 if there are no
 *           more sequences to read in the file.
 *           
 *           CloseFASTA() "always succeeds" and returns void.
 */
FASTAFILE *
OpenFASTA(char *seqfile)
{
  FASTAFILE *ffp;

  ffp = malloc(sizeof(FASTAFILE));
  ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
  if (ffp->fp == NULL) { free(ffp); return NULL; } 
  if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL)
    { free(ffp); return NULL; }
  return ffp;
}

int
ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, int *ret_L)
{
  char *s;
  char *name;
  char *seq;
  int   n;
  int   nalloc;
  
  /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
   */
  if (ffp->buffer[0] != '>') return 0;    

  /* Parse out the name: the first non-newline token after the >
   */
  s  = strtok(ffp->buffer+1, "\n");
  name = malloc(sizeof(char) * (strlen(s)+1));
  strcpy(name, s);

  /* Everything else 'til the next descline is the sequence.
   * Note the idiom for dynamic reallocation of seq as we
   * read more characters, so we don't have to assume a maximum
   * sequence length.
   */
  seq = malloc(sizeof(char) * 128);     /* allocate seq in blocks of 128 residues */
  nalloc = 128;
  n = 0;
  while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp))
    {
      if (ffp->buffer[0] == '>') break;	/* a-ha, we've reached the next descline */

      for (s = ffp->buffer; *s != '\0'; s++)
	{
	  if (! isalpha(*s)) continue;  /* accept any alphabetic character */

	  seq[n] = *s;                  /* store the character, bump length n */
	  n++;
	  if (nalloc == n)	        /* are we out of room in seq? if so, expand */
	    {			        /* (remember, need space for the final '\0')*/
	      nalloc += 128;
	      seq = realloc(seq, sizeof(char) * nalloc);
	    }
	}
    }
  seq[n] = '\0';

  *ret_name = name;
  *ret_seq  = seq;
  *ret_L    = n;
  return 1;
}      

void
CloseFASTA(FASTAFILE *ffp)
{
  fclose(ffp->fp);
  free(ffp);
}




/* what follows is a useful idiom: when you're writing a .c file that's supposed
 * to be a module of library functions, include one or more "test drivers".
 * These are small main()'s, normally ifdef'ed out of the code, that
 * enable the .c file to be compiled into one or more standalone test programs.
 * This lets you test your module in relative isolation, which tends
 * to lead to faster debugging and more robust code. It also
 * provides a convenient way to document a working minimal API: for example,
 * the main() here is a minimal FASTA reader. And it also tends to
 * have a useful psychological effect on you: it tends to encourage you
 * to simplify your APIs, so that small test programs can demonstrate
 * the full power of the API.
 */

#ifdef TEST_FASTA_STUFF
/* Test the fasta parsing API.
 *  to compile:  gcc -o test -DTEST_FASTA_STUFF -Wall -g fasta.c 
 *  to run:      ./test myseqs.fa 
 */
int
main(int argc, char **argv)
{
  FASTAFILE *ffp;
  char *seq;
  char *name;
  int   L;
				/* argv[1] is the name of a FASTA file */
  ffp = OpenFASTA(argv[1]);
  while (ReadFASTA(ffp, &seq, &name, &L))
    {
      printf(">%s\n", name);
      printf("%s\n",  seq);

      free(seq);
      free(name);
    }
  CloseFASTA(ffp);
  exit(0);
}
#endif /*TEST_FASTA_STUFF*/
