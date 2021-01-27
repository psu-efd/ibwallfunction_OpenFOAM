forAll (mesh.boundary(), patchI)
{
    if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
    {
        const immersedBoundaryFvPatch& ibPatch =
            refCast<const immersedBoundaryFvPatch>
            (
                mesh.boundary()[patchI]
            );
        immersedBoundaryFvPatch& ref_ibPatch =
            const_cast<immersedBoundaryFvPatch&>(ibPatch);

		ref_ibPatch.updateTriMesh();
	}	
}
