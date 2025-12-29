/*  MOLEKEL, Version 4.3, Date: 11.Nov.02
 *  Copyright (C) 2000-2002 Stefan Portmann (CSCS/ETHZ)
 *  (original IRIX GL implementation, concept and data structure
 *   by Peter F. Fluekiger, CSCS/UNI Geneva)
 *
 *  This software makes use of the 
 *  GLUT http://reality.sgi.com/mjk/glut3/glut3.html
 *  GLUI http://www.cs.unc.edu/~rademach/glui/
 *  libtiff http://www.libtiff.org/tiff-v3.5.5.tar.gz
 *  libjpeg ftp://ftp.uu.net/graphics/jpeg
 *  and in some versions of the 
 *  Mesa http://www.mesa3d.org/
 *  and the 
 *  libimage https://toolbox.sgi.com/toolbox/src/haeberli/libimage/index.html#dl
 *  libraries.
 *  An adapted version of the tr library by Brian Paul
 *  (http://www.mesa3d.org/brianp/TR.html)
 *  is part of the distribution.
 *
 *  The binary code is available free of charge but is not in the
 *  public domain. See license for details on conditions and restrictions.
 *  Source code is only available in the framework of a collaboration. 
 *
 *  Info: http://www.cscs.ch/molekel/
**/


char warranty_text[] = {"\nNO WARRANTY\n\
 \n\
3. DISCLAIMER OF WARRANTY. UNLESS SPECIFIED IN THIS\n\
AGREEMENT, ALL EXPRESS OR IMPLIED CONDITIONS,\n\
REPRESENTATIONS AND WARRANTIES, INCLUDING ANY IMPLIED\n\
WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR\n\
PURPOSE OR NON-INFRINGEMENT ARE DISCLAIMED, EXCEPT TO THE\n\
EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY\n\
INVALID. \n\
\n\
4. LIMITATION OF LIABILITY. TO THE EXTENT NOT PROHIBITED BY\n\
LAW, IN NO EVENT WILL CSCS OR ITS LICENSORS BE LIABLE FOR ANY\n\
LOST REVENUE, PROFIT OR DATA, OR FOR SPECIAL, INDIRECT,\n\
CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER\n\
CAUSED REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF\n\
OR RELATED TO THE USE OF OR INABILITY TO USE SOFTWARE, EVEN IF\n\
CSCS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES. In no\n\
event will CSCS's liability to you, whether in contract, tort \n\
(including negligence), or otherwise, exceed the amount paid by \n\
you for Software under this Agreement. The foregoing limitations \n\
will apply even if the above stated warranty fails of its \n\
essential purpose. \n\
\n"};

char license_text[] = {"\n\
 *  This software makes use of the\n\
 *  GLUT http://reality.sgi.com/mjk/glut3/glut3.html\n\
 *  GLUI http://www.cs.unc.edu/~rademach/glui/\n\
 *  libtiff http://www.libtiff.org/tiff-v3.5.5.tar.gz\n\
 *  and in some versions of the \n\
 *  Mesa http://www.mesa3d.org/\n\
 *  and the \n\
 *  libimage https://toolbox.sgi.com/toolbox/src/haeberli/libimage/index.html#dl\n\
 *  libraries.\n\
\n\
CSCS Binary Code License Agreement\n\
----------------------------------\n\
\n\
1. LICENSE TO USE. CSCS grants you a non-exclusive and\n\
non-transferable license for the internal use only of\n\
the accompanying software and documentation and\n\
any error corrections provided by CSCS (collectively \"Software\"). \n\
\n\
2. RESTRICTIONS Software is confidential and copyrighted.\n\
Title to Software and all associated INTELLECTUAL PROPERTY RIGHTS\n\
is retained by CSCS and/or its licensors. Except as specifically\n\
authorized by any Official Document from CSCS, you may not make\n\
copies of Software, other than a single copy of Software for\n\
archival purposes. Unless enforcement is prohibited by applicable\n\
law, you may not modify, decompile. No right, title or interest\n\
in or to any trademark, service mark, logo or trade name of CSCS\n\
or its licensors is granted under this Agreement. \n\
\n\
3. DISCLAIMER OF WARRANTY. UNLESS SPECIFIED IN THIS\n\
AGREEMENT, ALL EXPRESS OR IMPLIED CONDITIONS,\n\
REPRESENTATIONS AND WARRANTIES, INCLUDING ANY IMPLIED\n\
WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR\n\
PURPOSE OR NON-INFRINGEMENT ARE DISCLAIMED, EXCEPT TO THE\n\
EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY\n\
INVALID. \n\
\n\
4. LIMITATION OF LIABILITY. TO THE EXTENT NOT PROHIBITED BY\n\
LAW, IN NO EVENT WILL CSCS OR ITS LICENSORS BE LIABLE FOR ANY\n\
LOST REVENUE, PROFIT OR DATA, OR FOR SPECIAL, INDIRECT,\n\
CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER\n\
CAUSED REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF\n\
OR RELATED TO THE USE OF OR INABILITY TO USE SOFTWARE, EVEN IF\n\
CSCS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.\n\
In no event will CSCS's liability to you, whether in contract,\n\
tort (including negligence), or otherwise, exceed the amount paid\n\
by you for Software under this Agreement. The foregoing limitations\n\
will apply even if the above stated warranty fails of its essential\n\
purpose. \n\
\n\
5. Termination. This Agreement is effective until terminated.\n\
You may terminate this Agreement at any time by destroying all\n\
copies of Software. This Agreement will terminate immediately\n\
without notice from CSCS if you fail to comply with any\n\
provision of this Agreement. Upon Termination, you must destroy\n\
all copies of Software. \n\
\n\
6. Severability. If any provision of this Agreement is held to be\n\
unenforceable, this Agreement will remain in effect with the\n\
provision omitted, unless omission would frustrate the intent of the\n\
parties, in which case this Agreement will immediately terminate. \n\
\n\
7. Integration. This Agreement is the entire agreement between you\n\
and CSCS relating to its subject matter. It supersedes all prior or\n\
contemporaneous oral or written communications, proposals,\n\
representations and warranties and prevails over any conflicting or\n\
additional terms of any quote, order, acknowledgment, or other\n\
communication between the parties relating to its subject matter\n\
during the term of this Agreement. No modification of this Agreement\n\
will be binding, unless in writing and signed by an authorized\n\
representative of each party. \n\
\n\
For inquiries please contact:\n\
Centro Svizzero di Calcolo Scientifico (CSCS)\n\
Swiss Federal Institute of Technology\n\
Via Cantonale, Galleria 2\n\
CH-6928 Manno\n\
Switzerland\n\
\n\
\n\
Info: http://www.cscs.ch/molekel/\n"};

