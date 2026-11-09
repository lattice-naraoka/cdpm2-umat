! Standalone material-point tests ONLY. Not an Abaqus utility replacement.
subroutine MutexInit(id)
 integer::id
end subroutine
subroutine MutexLock(id)
 integer::id
end subroutine
subroutine MutexUnlock(id)
 integer::id
end subroutine
subroutine GETRANK(rank)
 integer::rank
 rank=0
end subroutine
subroutine GETNUMCPUS(n)
 integer::n
 n=1
end subroutine
subroutine GETPARTINFO(id,kind,name,localid,rc)
 integer::id,kind,localid,rc
 character(*)::name
 name='TEST';localid=id;rc=0
end subroutine
subroutine GETINTERNAL(name,id,kind,internalid,rc)
 integer::id,kind,internalid,rc
 character(*)::name
 internalid=id;rc=0
end subroutine
subroutine xit()
 stop 99
end subroutine
