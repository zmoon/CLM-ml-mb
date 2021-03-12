program CLMml

  use decompMod, only : bounds_type, get_clump_bounds
  use ColumnType, only : col
  use PatchType, only : patch
  use clm_instMod, only : clm_instInit
  use filterMod, only : filter, allocFilters
  use CLMml_driver, only : CLMml_drv

  type(bounds_type) :: bounds

  ! Define grid cell (g), land unit (l), column (c), and patch (p) bounds
  ! for CLM g/l/c/p hierarchy. CLM processes clumps of gridcells (and
  ! associated subgrid-scale entities) with length defined by
  ! begg/endg, begl/endl, begc/endc, and begp/endp. This code assumes
  ! that a grid cell has one land unit with one column and one patch. It
  ! processes a single grid cell.

  print *, "Getting CLM g/l/c/p ..."
  call get_clump_bounds (bounds)

  ! Initialize instances of all derived types

  print *, "Initializing derived types ..."
  call col%Init (bounds%begc, bounds%endc)
  call patch%Init (bounds%begp, bounds%endp)
  call clm_instInit (bounds)

  ! Allocate filters

  print *, "Allocating filters ..."
  call allocFilters (filter, bounds%begp, bounds%endp, bounds%begc, bounds%endc)

  ! Run model

  print *, "Running the model ..."
  call CLMml_drv (bounds)

end program CLMml
