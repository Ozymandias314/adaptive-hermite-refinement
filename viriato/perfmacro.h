  !. FPP macro to simplify RZG PERFLIB calls (performance-counter)
# ifdef perf
#    define Myperfinit   call Perfinit
#    define Myperfon(X)  call Perfon(X)
#    define Myperfout(X) call Perfout(X)
#    define Myperfoff    call Perfoff
# else
#    define Myperfinit
#    define Myperfon(X)
#    define Myperfoff
#    define Myperfout(X)
# endif
