!!!!!!!!
!
!fichier g2s_cast.md9
!
!  module cast
!  ###########
! permet d'interpreter des entiers en chaines de caracteres et vice-versa
!  - c2i : transforme un caractere en chiffre
!  - s2i : transforme une chaine de caracteres en nombre
!  - i2s : transforme un nombre en chaine de caracteres
!
!
!!!!!!!!


module cast
  implicit none

  integer,parameter :: izero=ichar('0')

contains

  function c2i(c)
    character,intent(in) :: c
    integer              :: c2i

    c2i = ichar(c) - izero
    if (c2i<0 .or. c2i>9) &
         & call XABORT("G2S: internal error, bad cast from string to integer")
  end function c2i

  function s2i(s)
    character(len=*),intent(in) :: s
    integer                     :: s2i
    integer   :: i,l,n

    s2i = 0
    l = len_trim(s)
    do i = 1,l
       n = c2i(s(i:i))
       s2i = s2i + n*10**(l-i)
    end do
  end function s2i

  function i2s(i)
    integer,intent(in) :: i
    character*12       :: i2s
    
    character*12 :: tmp
    
    tmp = ' '
    write(tmp,'(i12)') i
    i2s = adjustl(tmp)
  end function i2s

end module cast
