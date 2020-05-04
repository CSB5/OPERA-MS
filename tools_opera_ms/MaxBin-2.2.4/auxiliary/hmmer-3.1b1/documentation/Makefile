top_srcdir     = ..
srcdir         = .


QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR $$subdir; \
		        ${MAKE} -C $$subdir
endif


.PHONY: pdf install uninstall clean distclean

pdf:
	${QUIET_SUBDIR0}userguide ${QUIET_SUBDIR1} pdf

install: 
	${QUIET_SUBDIR0}man       ${QUIET_SUBDIR1} install

uninstall:
	${QUIET_SUBDIR0}man       ${QUIET_SUBDIR1} uninstall

clean:
	${QUIET_SUBDIR0}man       ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}userguide ${QUIET_SUBDIR1} clean
	-rm -f *~

distclean:
	${QUIET_SUBDIR0}man       ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}userguide ${QUIET_SUBDIR1} distclean
	-rm -f Makefile

################################################################
# HMMER - Biological sequence analysis with profile HMMs
# Version 3.1b1; May 2013
# Copyright (C) 2013 Howard Hughes Medical Institute.
# Other copyrights also apply. See the COPYRIGHT file for a full list.
# 
# HMMER is distributed under the terms of the GNU General Public License
# (GPLv3). See the LICENSE file for details.
# 
# SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/branches/3.1/documentation/Makefile.in $
# SVN $Id: Makefile.in 3649 2011-08-11 15:12:30Z eddys $
################################################################
