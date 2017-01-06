module types
  
  implicit none

  type typemod
     real,    allocatable :: vp(:, :), vs(:, :), ro(:, :)
     real,    allocatable :: vpe(:, :), vse(:, :), roe(:, :)
     real,    allocatable :: bux(:, :), buz(:, :)
     real,    allocatable :: mu0(:, :), mue(:, :)
     real,    allocatable :: lb0(:, :), lbmu(:, :)
  end type typemod

end module types
