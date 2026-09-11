module core
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, core!"
  end subroutine say_hello
end module core
