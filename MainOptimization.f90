program MainOptimization
    ! 1. Librerias usadas
    use Optimization_module
    use Paraview_Module
    implicit none
    type (Optimization)                                  :: OptimizationModel
    type(PostprocessingInfo)                             :: ParaviewModel
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                          TOPOLOGICAL OPTIMIZATION CODE USING FEM IN LINEAR ELASTICITY                      !
    !                                                   (TOP-L)                                                  !
    !                                                                                                            !
    ! Note: This software is open source and is provided for academic teaching and research purposes only. It    !
    !       may be freely used, modified and distributed. No warranty or technical support is provided for it.   !
    !       Users are expected to understand and take full responsibility for its use. If used for research      !
    !       purposes, please cite the repository.                                                                !
    !                                                                                                            !
    ! Important: The TOP-L program uses the HSL-LA87 library, so it is necessary to have the LaPack and BLAS     !
    !       libraries (preferably the optimized versions) installed and the Metis library.                       !
    !                                                                                                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2. Ingreso de parametros y propiedades de estructura
    write(unit=*, fmt=*) '1. Structure Properties'
    ! 2.1. Mechanical parameters
    call SetAnalysisType(OptimizationModel,'SolidIso')               ! Set PlaneStress(2D), PlainStrain(2D), SolidIso(3D)
    call SetElementType(OptimizationModel,'hexa8')                   ! tria3 tria6, quad4, quad8, tetra4, tetra10, hexa8, hexa20
    call SetThickness(OptimizationModel,0.0d0)                       ! Only for 2D analysis
    call SetYoungModulus(OptimizationModel,210000.0d0)                 ! Young modulus
    call SetPoissonModulus(OptimizationModel,0.3d0)                  ! Poisson modulus
    call SetGaussAprox(OptimizationModel,3)                          ! Can use up to 5 gauss points
    ! 2.2. Optimization parameters
    call SetOptimizationAlgorithm(OptimizationModel,'OCM')           ! Optimality criteria (OCM), Method of moving Asym. (MMA)
    call SetMinDensityValue(OptimizationModel,1.0d-2)                ! Min. Density value
    call SetVolFraction(OptimizationModel,0.5d0)                     ! Limiting volume fraction (lower limit of TOP)
    call SetMutationRate(OptimizationModel,0.2d0)                    ! Rate of change/mutation of optimization
    call SetFilterRadius(OptimizationModel,10.0d0*1.5)               ! Smoothing filter radius (1.5 times de FE size)
    call SetMaxIterations(OptimizationModel,30)                      ! Maximum number of iterations
    call SetPenalFactor(OptimizationModel,3.0d0)                     ! SIMP method penalty factor
    call SetOptimizationTolerance(OptimizationModel,0.01d0)            ! Optimization tolerance
    call ReadFiles(OptimizationModel)
    ! 3. Applying the topology optimization
    write(unit=*, fmt=*) '2. Topology optimization'
    call TopologyOptimizationProcess(OptimizationModel)
    ! 4. Postprocesado con paraview
    write(unit=*, fmt=*) '3. PostProcessing'
    call SetPostProcesingFilter(OptimizationModel,0.5d0)               ! PostProcessing filter (to eliminate inter. densities)
    ! note: leave whatever you want to postprocess and generate the 
    !       file for viewing with paraview. If you don't want to 
    !       postprocess just comment the line.
    ! Displacement, Strain, Stress, StrainEnergy, DensityDist
    call PlottingSelection(ParaviewModel,'DensityDist')
    call GetParaviewFiles(OptimizationModel,ParaviewModel)
end program MainOptimization