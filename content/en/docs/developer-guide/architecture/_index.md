---
title: "Architecture"
linkTitle: "Architecture"
weight: 1
description: >
  Overview of the multiPhysicsFoam modular architecture.
---

The class diagram below shows the core classes of multiPhysicsFoam and how they relate to each other and to the OpenFOAM base classes.

## UML Class Diagram

```mermaid
classDiagram
    %% =========================================================
    %% OpenFOAM Base Classes (simplified stubs)
    %% =========================================================
    class IOdictionary {
        <<OpenFOAM>>
    }
    class MeshObject_fvMesh {
        <<OpenFOAM>>
        +movePoints() bool
        +updateMesh(mapPolyMesh) void
    }
    class Pair_word {
        <<OpenFOAM>>
    }
    class IOobject {
        <<OpenFOAM>>
    }
    class PtrList_T {
        <<OpenFOAM>>
    }
    class fixedValueFvPatchField_T {
        <<OpenFOAM>>
        +updateCoeffs() void
        +evaluate() void
    }
    class fixedGradientFvPatchField_T {
        <<OpenFOAM>>
        +updateCoeffs() void
        +evaluate() void
    }
    class coupledFvPatchField_T {
        <<OpenFOAM>>
    }
    class ggiLduInterfaceField {
        <<foam-extend>>
        +initInterfaceMatrixUpdate() void
        +updateInterfaceMatrix() void
    }

    %% =========================================================
    %% interfaceKey
    %% =========================================================
    class interfaceKey {
        -bool ordered_
        +bool ordered() const
        +hash operator(interfaceKey) label
    }
    Pair_word <|-- interfaceKey

    %% =========================================================
    %% globalPolyPatch
    %% =========================================================
    class globalPolyPatch {
        -polyMesh& mesh_
        -word patchName_
        -polyPatch& patch_
        -standAlonePatch* globalPatchPtr_
        -labelList* pointToGlobalAddrPtr_
        -labelList* faceToGlobalAddrPtr_
        -labelList* globalMasterToCurrentProcPointAddrPtr_
        -PrimitivePatchInterpolation* interpPtr_
        +word patchName() const
        +standAlonePatch& globalPatch() const
        +PrimitivePatchInterpolation& interpolator() const
        +labelList& pointToGlobalAddr() const
        +labelList& faceToGlobalAddr() const
        +labelList& globalMasterToCurrentProcPointAddr() const
        +Field_T patchPointToGlobal(Field_T)
        +Field_T globalPointToPatch(Field_T)
        +Field_T patchFaceToGlobal(Field_T)
        +Field_T globalFaceToPatch(Field_T)
    }

    %% =========================================================
    %% interfaceToInterfaceMapping  (abstract, runtime-selectable)
    %% =========================================================
    class interfaceToInterfaceMapping {
        <<abstract>>
        #dictionary dict_
        #primitivePatch& patchA_
        #primitivePatch& patchB_
        #globalPolyPatch& globalPatchA_
        #globalPolyPatch& globalPatchB_
        +transferFacesZoneToZone_scalar()* void
        +transferFacesZoneToZone_vector()* void
        +transferPointsZoneToZone_scalar()* void
        +transferPointsZoneToZone_vector()* void
        +New()$ autoPtr_interfaceToInterfaceMapping
    }
    class amiInterfaceToInterfaceMapping {
        <<AMI>>
        -autoPtr_amiZoneInterpolation interpolatorPtr_
        -List_labelPair* zoneAPointAddressingPtr_
        -FieldField_scalar* zoneAPointWeightsPtr_
        -List_labelPair* zoneBPointAddressingPtr_
        -FieldField_scalar* zoneBPointWeightsPtr_
        +transferFacesZoneToZone_T() void
        +transferPointsZoneToZone_T() void
    }
    class directMapInterfaceToInterfaceMapping {
        <<directMap>>
        -autoPtr_labelIOList zoneAToZoneBFaceMapPtr_
        -autoPtr_labelIOList zoneBToZoneAFaceMapPtr_
        -autoPtr_labelIOList zoneAToZoneBPointMapPtr_
        -autoPtr_labelIOList zoneBToZoneAPointMapPtr_
        -autoPtr_scalarField minEdgeLengthsPtr_
        +transferFacesZoneToZone_T() void
        +transferPointsZoneToZone_T() void
    }
    class rbfInterfaceToInterfaceMapping {
        <<RBF, needs Eigen>>
        -shared_ptr_RBFInterpolation zoneAToZoneBInterpolatorPtr_
        -shared_ptr_RBFInterpolation zoneBToZoneAInterpolatorPtr_
        +transferFacesZoneToZone_T() void
        +transferPointsZoneToZone_T() void
    }
    class amiZoneInterpolation {
        <<helper>>
    }
    interfaceToInterfaceMapping <|-- amiInterfaceToInterfaceMapping
    interfaceToInterfaceMapping <|-- directMapInterfaceToInterfaceMapping
    interfaceToInterfaceMapping <|-- rbfInterfaceToInterfaceMapping
    amiInterfaceToInterfaceMapping *-- amiZoneInterpolation : interpolatorPtr_
    interfaceToInterfaceMapping "1" o-- "2" globalPolyPatch : globalPatchA_ / globalPatchB_

    %% =========================================================
    %% regionInterfaceType  (abstract base for all interfaces)
    %% =========================================================
    class regionInterfaceType {
        <<abstract>>
        #IOdictionary multiRegionProperties_
        #IOdictionary regionInterfaceProperties_
        #Time& runTime_
        #fvPatch& patchA_
        #fvPatch& patchB_
        #fvMesh& meshA_
        #fvMesh& meshB_
        #autoPtr_globalPolyPatch globalPatchAPtr_
        #autoPtr_globalPolyPatch globalPatchBPtr_
        #Switch attachedA_ / attachedB_
        #bool changing_
        #bool moving_
        #Switch coupled_
        #scalar couplingStartTime_
        #label interpolatorUpdateFrequency_
        #autoPtr_faMesh aMeshPtr_
        #wordList curvatureCorrectedSurfacePatches_
        #autoPtr_interfaceToInterfaceMapping interfaceToInterfacePtr_
        +correct()* void
        +getMinDeltaT()* scalar
        +info()* void
        +interfaceName() word
        +updateInterpolatorAndGlobalPatches() void
        +attach() void
        +detach() void
        +update() void
        +updateK() void
        +movePoints() bool
        +updateMesh(mapPolyMesh) void
        +couplesField(word) bool
        +usesOneInterfaceRes() Switch
        +rawJumpRes_T(word) scalarField
        +rawFluxRes_T(word) scalarField
        +rawIntRes(word) scalarField
        +normJumpRes_T(word) scalar
        +normFluxRes_T(word) scalar
        +normIntRes(word) scalar
        +interpolateFacesFromA_T(Field_T) Field_T
        +interpolateFacesFromB_T(Field_T) Field_T
        +interpolatePointsFromA_T(Field_T) Field_T
        +interpolatePointsFromB_T(Field_T) Field_T
        +Switch& coupled()
        +faMesh& aMesh()
        +fvPatch& patchA() / patchB()
        +fvMesh& meshA() / meshB()
        +globalPolyPatch& globalPatchA() / globalPatchB()
        +interfaceToInterfaceMapping& interfaceToInterface()
    }
    IOdictionary           <|-- regionInterfaceType
    MeshObject_fvMesh      <|-- regionInterfaceType
    interfaceKey           <|-- regionInterfaceType
    regionInterfaceType "1" *-- "2" globalPolyPatch       : globalPatchAPtr_ / globalPatchBPtr_
    regionInterfaceType  *--  interfaceToInterfaceMapping : interfaceToInterfacePtr_

    %% =========================================================
    %% Concrete Interface Types
    %% =========================================================
    class heatTransferInterface {
        -dictionary dict_
        -dimensionedScalar jumpCoefficient_
        +jumpCoefficient() dimensionedScalar
        +correct() void
        +getMinDeltaT() scalar
        +info() void
    }
    class capillaryInterface {
        -dictionary dict_
        -dimensionedScalar sigma0_
        -areaScalarField sigma_
        -autoPtr_areaVectorField UsPtr_
        -autoPtr_edgeScalarField phisPtr_
        +Us() areaVectorField&
        +Phis() edgeScalarField&
        +sigma() areaScalarField&
        +surfaceTensionForce() vectorField
        +tangentialSurfaceTensionForce() vectorField
        +normalSurfaceTensionForce() vectorField
        +totalSurfaceTensionForce() vector
        +totalViscousForce() vector
        +totalPressureForce() vector
        +updateUs() void
        +updatePhis() void
        +correct() void
        +getMinDeltaT() scalar
        +info() void
    }
    class massTransferInterface {
        -dictionary dict_
        -dimensionedScalar hlv0_
        -areaScalarField hlv_
        -dimensionedScalar TSat0_
        -dimensionedScalar mDots0_
        -areaScalarField mDotInterface
        -autoPtr_areaScalarField mDotsPtr_
        -autoPtr_areaScalarField fluxMTJumpPtr_
        +mDotS() areaScalarField&
        +fluxMTJump() areaScalarField&
        +hlv() areaScalarField&
        +saturatedFlux() scalarField
        +updateMDotS() void
        +correct() void
        +getMinDeltaT() scalar
        +info() void
    }
    class fsiInterface {
        -dictionary dict_
        -vector totalForce_
        +makeInterfaceToInterface() void
        +setTotalForce(vector) void
        +rawIntRes(word) scalarField
        +normIntRes(word) scalar
        +correct() void
        +getMinDeltaT() scalar
        +info() void
    }
    class defaultInterface {
        +correct() void
        +getMinDeltaT() scalar
        +info() void
    }
    regionInterfaceType <|-- heatTransferInterface
    regionInterfaceType <|-- capillaryInterface
    regionInterfaceType <|-- massTransferInterface
    regionInterfaceType <|-- fsiInterface
    regionInterfaceType <|-- defaultInterface

    %% =========================================================
    %% interfaceToInterfaceCoupleManager
    %% (mixin base for coupled boundary conditions)
    %% =========================================================
    class interfaceToInterfaceCoupleManager {
        #fvPatch& patch_
        #fvMesh& localRegion_
        #wordList typeNames_
        #word neighbourFieldName_
        #rgInterface(word) regionInterfaceType&
        #nbrMesh() fvMesh&
        #nbrPatch() fvPatch&
        #interfaceCoupled() Switch
        +interpolateFromNbrField_T(Field_T) Field_T
        +interpolateToNbrField_T(Field_T) Field_T
        +refMesh() fvMesh&
        +refPatch() fvPatch&
        +refPatchID() label
        +neighbourFieldName() word
        +neighbourRegion() fvMesh&
        +neighbourPatchID() label
        +neighbourPatch() fvPatch&
        +neighbourPatchField_T() fvPatchField_T&
        +updateRegionInterface() void
    }
    interfaceToInterfaceCoupleManager --> regionInterfaceType : rgInterface() looks up\nfrom mesh object registry

    %% =========================================================
    %% Boundary Conditions — Partitioned (DNA) coupling
    %% =========================================================
    class genericRegionCoupledFluxFvPatchField_T {
        <<Partitioned / fixedGradient>>
        -word diffusivityName_
        -autoPtr_accelerationModel_T accModel_
        -Switch nonOrthCorr_
        -Switch secondOrder_
        #useDirectFlux() bool
        +fluxJump() Field_T
        +flux() Field_T
        +diffusivityName() word
        +updateCoeffs() void
        +evaluate() void
        +rawResidual() scalarField
        +normResidual() scalar
        +ofNormResidual() scalar
        +diffCoeff() scalarField
    }
    class genericRegionCoupledJumpFvPatchField_T {
        <<Partitioned / fixedValue>>
        -word diffusivityName_
        -word jumpCoefficientName_
        -autoPtr_accelerationModel_T accModel_
        -Switch nonOrthCorr_
        -Switch secondOrder_
        #useDirectValue() bool
        #updatePhi() void
        +valueJump() Field_T
        +value() Field_T
        +diffusivityName() word
        +updateCoeffs() void
        +flux() Field_T
        +snGrad() Field_T
        +gradientBoundaryCoeffs() Field_T
        +rawResidual() scalarField
        +normResidual() scalar
        +ofNormResidual() scalar
        +diffCoeff() scalarField
    }
    fixedGradientFvPatchField_T        <|-- genericRegionCoupledFluxFvPatchField_T
    interfaceToInterfaceCoupleManager  <|-- genericRegionCoupledFluxFvPatchField_T
    fixedValueFvPatchField_T           <|-- genericRegionCoupledJumpFvPatchField_T
    interfaceToInterfaceCoupleManager  <|-- genericRegionCoupledJumpFvPatchField_T
    genericRegionCoupledFluxFvPatchField_T *-- accelerationModel_T : accModel_
    genericRegionCoupledJumpFvPatchField_T *-- accelerationModel_T : accModel_

    %% concrete partitioned BCs
    class regionCoupledFluxBCs {
        <<group: genericRegionCoupledFlux>>
        regionCoupledScalarFlux
        regionCoupledHeatFlux
        regionCoupledPhaseChangeHeatFlux
        regionCoupledPressureFlux
        regionCoupledVelocityFlux
        regionCoupledMassTransferVelocityFlux
        regionCoupledTraction
    }
    class regionCoupledJumpBCs {
        <<group: genericRegionCoupledJump>>
        regionCoupledScalarJump
        regionCoupledTemperatureJump
        regionCoupledPressureValue
        regionCoupledMassTransferPressureValue
        regionCoupledVelocityValue
        regionCoupledMassTransferVelocityValue
    }
    genericRegionCoupledFluxFvPatchField_T <|-- regionCoupledFluxBCs
    genericRegionCoupledJumpFvPatchField_T <|-- regionCoupledJumpBCs

    %% =========================================================
    %% Convergence acceleration (runtime-selectable)
    %% =========================================================
    class accelerationModel_T {
        <<abstract>>
        #Time& runTime_
        #scalar curTime_
        #label corr_
        #Field_T prevFld_
        #Field_T resFld_
        #Field_T prevResFld_
        #scalar initRelax_
        +initialize(Field_T) void
        +relax(Field_T)* void
        +New()$ autoPtr_accelerationModel
    }
    class fixedRelaxation { <<fixed>> }
    class aitkenRelaxation { <<aitken>> }
    class IQNILS { <<IQNILS>> }
    accelerationModel_T <|-- fixedRelaxation
    accelerationModel_T <|-- aitkenRelaxation
    accelerationModel_T <|-- IQNILS

    %% =========================================================
    %% Boundary Conditions — Monolithic coupling
    %% NOT COMPILED: commented out in Make/files ("NEED TO BE PORTED
    %% STILL"), still based on foam-extend's ggiLduInterfaceField
    %% =========================================================
    class monolithicCouplingFvPatchField_T {
        <<Monolithic / not ported>>
        -regionCoupleFvPatch& regionCouplePatch_
        -word remoteFieldName_
        -scalarField matrixUpdateBuffer_
        -Field_T originalPatchField_
        -label curTimeIndex_
        +remoteFieldName() word
        +setRemoteFieldName(word) void
        +setMatrixUpdateBuffer(scalarField) void
        +shadowPatchField() monolithicCouplingFvPatchField_T&
        +patchNeighbourField() Field_T
        +updateCoeffs() void
        +snGrad() Field_T
        +initInterfaceMatrixUpdate() void
        +updateInterfaceMatrix() void
        +transformCoupleField(scalarField, direction) void
        +doTransform() bool
    }
    ggiLduInterfaceField              <|-- monolithicCouplingFvPatchField_T
    coupledFvPatchField_T             <|-- monolithicCouplingFvPatchField_T
    interfaceToInterfaceCoupleManager <|-- monolithicCouplingFvPatchField_T

    %% =========================================================
    %% regionType  (abstract base for all physics regions)
    %% =========================================================
    class regionType {
        <<abstract>>
        #refPtr_dynamicFvMesh meshPtr_
        #HashTable_scalar maxCorr_
        #HashTable_scalar relativeTolerance_
        #HashTable_scalar convergenceTolerance_
        #HashTable_int convergenceInfoFeq_
        #HashPtrTable_fvMatrix_scalar fvScalarMatrices
        #HashPtrTable_fvMatrix_vector fvVectorMatrices
        #HashPtrTable_fvMatrix_symmTensor fvSymmTensorMatrices
        #HashPtrTable_fvMatrix_tensor fvTensorMatrices
        #lookupOrRead_T(fvMesh, word, ...) refPtr_T
        +regionTypeName()* word
        +usesPIMPLE()* bool
        +correct()* void
        +solveRegion()* void
        +prePredictor()* void
        +momentumPredictor()* void
        +pressureCorrector()* void
        +meshMotionCorrector()* void
        +setCoupledEqns(word)* void
        +selfCoupledFieldNames() wordList
        +postSolve()* void
        +getMinDeltaT()* scalar
        +regionName() word
        +mesh() dynamicFvMesh&
        +runTime() Time&
        +relaxEqn_T(fvMatrix_T) void
        +getCoupledEqn_T(word) fvMatrix_T&
        +clearCoupledEqn_T(T, word) bool
        +foundCoupledEqn(word) bool
        +getObject_T(word, fvMesh) T&
        +usesPicard(word) bool
        +converged(word, int, scalar, scalar) bool
    }
    IOdictionary <|-- regionType

    %% concrete region types compiled in Make/files
    %% (buoyantPimpleFluid, pUCoupledIcoFluid, interTrackFluid and
    %% turbulentTransportTemperature exist in src/regions but are disabled)
    class fluidRegionTypes {
        <<group>>
        icoFluid
        pimpleFluid
        rhoPimpleFluid
        overPimpleFluid
        interFluid
        laplacianFluid
        potentialFluid
    }
    class thermalRegionTypes {
        <<group>>
        conductTemperature
        transportTemperature
        conductPotentialTemperature
        NTGK
    }
    class solidRegionTypes {
        <<group>>
        solidStVK
    }
    regionType <|-- fluidRegionTypes
    regionType <|-- thermalRegionTypes
    regionType <|-- solidRegionTypes

    %% =========================================================
    %% regionTypeList / regionInterfaceTypeList
    %% =========================================================
    class regionTypeList {
        -IOdictionary dict_
        -Time& runTime_
        -List_Tuple2 regionProperties_
        -wordList meshNames_
        +usesPIMPLE(bool) bool
        +selfCoupledFldNames() hashedWordList
        +reset(List_Tuple2) void
        +preSolve() void
        +solveRegion() void
        +solvePIMPLE() void
        +meshMotionCorrector() void
        +setCoupledEqns(word) void
        +postSolve() void
        +postSolvePIMPLE() void
        +getMinDeltaT() scalar
        +clear() void
    }
    PtrList_T <|-- regionTypeList : PtrList~regionType~
    regionTypeList "1" *-- "1..*" regionType : elements

    class regionInterfaceTypeList {
        -label index_
        -Time& runTime_
        -HashTable partitionedCoupledFields_
        -HashTable monolithicCoupledFields_
        -hashedWordList pcFldNames_
        -hashedWordList mcFldNames_
        -regionInterfaceProperties partitionedTypeInterfaces_
        -regionInterfaceProperties monolithicTypeInterfaces_
        +partitionedCoupledFields() HashTable
        +monolithicCoupledFields() HashTable
        +pcFldNames() hashedWordList
        +mcFldNames() hashedWordList
        +attach() void
        +detach() void
        +update() void
        +getMinDeltaT() scalar
        +reset(regionInterfaceProperties) void
        +setFieldNamesPartitionedCoupling(regionInterfaceProperties) void
        +setFieldNamesMonolithicCoupling(regionInterfaceProperties) void
    }
    PtrList_T <|-- regionInterfaceTypeList : PtrList~regionInterfaceType~
    regionInterfaceTypeList "1" *-- "1..*" regionInterfaceType : elements

    %% =========================================================
    %% dnaControl  (Dirichlet-Neumann Alternating convergence)
    %% =========================================================
    class dnaControl {
        -word fldName_
        -IOdictionary multiRegionProperties_
        -regionInterfaceTypeList& interfaces_
        -label maxCoupleIter_
        -label corr_
        -List_fieldData dnaResControl_
        -List_scalar globalMaxIntRes_
        -List_scalar globalMaxJumpRes_
        -List_scalar globalMaxFluxRes_
        -HashPtrTable_volScalarField initIntResFlds_ / initJumpResFlds_ / initFluxResFlds_
        -HashPtrTable_volScalarField finalIntResFlds_ / finalJumpResFlds_ / finalFluxResFlds_
        -createResFlds() void
        +read() void
        +loop() bool
        +criteriaSatisfied() bool
        +corr() label
        +dict() dictionary
    }
    IOobject <|-- dnaControl
    dnaControl "1" o-- "1" regionInterfaceTypeList : monitors residuals

    %% =========================================================
    %% multiRegionSystem  (top-level orchestrator)
    %% =========================================================
    class multiRegionSystem {
        -Time& runTime_
        -autoPtr_regionTypeList regions_
        -autoPtr_regionInterfaceTypeList interfaces_
        -hashedWordList partitionedCoupledFldNames_
        -HashPtrTable_dnaControl dnaControls_
        +preSolve() void
        +solve() void
        +postSolve() void
        +getMinDeltaT() scalar
        +regions() regionTypeList&
        +interfaces() regionInterfaceTypeList&
        -assembleCoupledFields_T(PtrList_volField_T, hashedWordList) void
        -assembleAndSolveEqns_M_T(word) void
    }
    IOdictionary         <|-- multiRegionSystem
    multiRegionSystem "1" *-- "1"   regionTypeList          : regions_
    multiRegionSystem "1" *-- "1"   regionInterfaceTypeList : interfaces_
    multiRegionSystem "1" *-- "1..*" dnaControl             : dnaControls_
```
