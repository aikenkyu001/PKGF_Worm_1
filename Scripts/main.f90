program pkgf_worm_main
    use pkgf_constants
    use pkgf_worm_engine_mod
    use pkgf_math_utils
    implicit none
    type(PKGFWormEngine) :: engine
    integer :: step, i, log_unit
    real(DP) :: px, py, avg_h, avg_t, avg_v
    
    print *, " === PKGF-Worm Sync Verification (Fortran) ==="
    call engine%init()
    
    ! ログディレクトリへ保存
    open(newunit=log_unit, file="../Logs/worm_log_fortran.txt", status="replace")
    write(log_unit, '(A)') "Step,HeadX,HeadY,Hunger,Tension,V,DetK"
    
    do step = 0, 300
        px = engine%nodes(1)%pos(1)
        py = engine%nodes(1)%pos(2)
        avg_h = 0.0d0; avg_t = 0.0d0; avg_v = 0.0d0
        do i = 1, N_NODES
            avg_h = avg_h + engine%nodes(i)%hunger()
            avg_t = avg_t + engine%nodes(i)%tension()
            avg_v = avg_v + sqrt(sum(engine%nodes(i)%v**2))
        end do
        avg_h = avg_h / real(N_NODES, DP)
        avg_t = avg_t / real(N_NODES, DP)
        avg_v = avg_v / real(N_NODES, DP)
        
        ! 厳密なCSV出力
        write(log_unit, '(I5,A,F18.10,A,F18.10,A,F18.10,A,F18.10,A,F18.10,A,F18.10)') &
            step, ",", px, ",", py, ",", avg_h, ",", avg_t, ",", avg_v, ",", determinant(engine%nodes(1)%K)
            
        if (mod(step, 50) == 0) then
            write(*, '(A,I5,A,F10.2,A,F10.2,A,F10.6,A,F10.6)') &
                "Step", step, " | Pos: (", px, ",", py, ") | Hunger:", avg_h, " | V:", avg_v
        end if
        
        call engine%step(step)
    end do
    
    close(log_unit)
    print *, " Log saved to Logs/worm_log_fortran.txt"
end program pkgf_worm_main
