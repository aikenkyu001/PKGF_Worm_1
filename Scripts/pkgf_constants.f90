module pkgf_constants
    implicit none
    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: N_NODES = 302
    integer, parameter :: DIM = 32
    real(DP), parameter :: PI = 3.14159265358979323846d0
    real(DP), parameter :: DT = 0.02d0
    ! 幾何学的定数: ALPHA は多様体の曲率結合定数として一元化
    real(DP), parameter :: ALPHA = 1.0d0 / real(DIM, DP) 
    ! 背景粘性: 物理世界の「重さ」
    real(DP), parameter :: VISCOSITY = 5.0d0 
end module pkgf_constants
