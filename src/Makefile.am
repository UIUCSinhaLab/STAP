AUTOMAKE_OPTIONS = foreign

#TODO : define values are still not passed from the configure script.
AM_CPPFLAGS = $(GSL_CFLAGS) -DStrongerSiteAttr=0 -DRandomGen=0

bin_PROGRAMS = seq2binding
seq2binding_SOURCES = seq2binding.cpp BindingPredictor.cpp BindingPredictor.h SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h
seq2binding_LDADD = $(GSL_LIBS)
