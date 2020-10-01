# -*- coding: utf-8 -*-
"""
Class for the material properties of a lamina
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

class Material():
    " An object for storing the lamina properties"

    def __init__(self,
                 E11=130e9, # Elastic modulus in the principal direction (Pa)
                 E22=9e9, # Elastic modulus in the transverse direction (Pa)
                 nu12=0.3, # Poisson's ratio
                 G12=4e9, # In-plane shear modulus (Pa)
                 density_volume=0, # Density (g/m3)
                 density_area=0,# Density (g/m2)
                 ply_t=0.0002): # Ply thickness (m)
        " Create an instance of the class Material"

        # Elastic modulus in the fibre direction (Pa)
        if not isinstance(E11, int)  and not isinstance(E11, float):
            raise MaterialDefinitionError(f"""
Attention, the elastic modulus must be a number! (E11 = {E11})""")
        if E11 <= 0:
            raise MaterialDefinitionError(f"""
Attention, the elastic modulus must be positive! (E11 = {E11})""")

        # Elastic modulus in the transverse direction (Pa)
        if not isinstance(E22, int) and not isinstance(E22, float):
            raise MaterialDefinitionError(f"""
Attention, the elastic modulus must be a number! (E22 = {E22})""")
        if E11 <= 0:
            raise MaterialDefinitionError(f"""
Attention, the elastic modulus must be positive! (E22 = {E22})""")

        # In-plane shear modulus (Pa)
        if not isinstance(G12, int) and not isinstance(G12, float):
            raise MaterialDefinitionError(f"""
Attention, the shear modulus must be a number! (G12 = {G12})""")
        if E11 <= 0:
            raise MaterialDefinitionError(f"""
Attention, the shear modulus must be positive! (G12 = {G12})""")

        # Poisson's ratio relating transverse deformation and axial loading
        if not isinstance(nu12, int) and not isinstance(nu12, float):
            raise MaterialDefinitionError(f"""
Attention, the Poisson\'s ratio must be a number! (nu12 = {nu12})""")

        # Density
        if not isinstance(density_volume, int) \
        and not isinstance(density_volume, float):
            raise MaterialDefinitionError(f"""
Attention, the density must be a number! (density_volume = {density_volume})""")
        if density_volume < 0:
            raise MaterialDefinitionError(f"""
Attention, the density must be positive! (density_volume = {density_volume})""")

        # Density
        if not isinstance(density_area, int) \
        and not isinstance(density_area, float):
            raise MaterialDefinitionError(f"""
Attention, the density must be a number! (density_area = {density_area})""")
        if density_area < 0:
            raise MaterialDefinitionError(f"""
Attention, the density must be positive! (density_area = {density_area})""")

        # Ply thickness (mm)
        if not isinstance(ply_t, int) \
        and not isinstance(ply_t, float):
            raise MaterialDefinitionError(f"""
Attention, the shear modulus must be a number! (ply_t = {ply_t})""")
        if ply_t <= 0:
            raise MaterialDefinitionError(f"""
Attention, the density must be positive! (ply_t= {ply_t})""")

        self.E11 = E11
        self.E22 = E22
        self.G12 = G12
        self.nu12 = nu12
        self.nu21 = nu12*E22/E11
        self.density_area = density_area
        self.density_volume = density_volume
        self.ply_t = ply_t

        if density_area == 0 and density_volume == 0:
            # give a value by default
            self.density_area = 300.5 # g/m^2
            self.density_volume = density_area / ply_t
        elif density_area == 0:
            self.density_area = density_volume * ply_t
        elif density_volume == 0:
            self.density_volume = density_area / ply_t

        # Lamina stiffnesses
        bla = 1 - nu12*self.nu21
        self.Q11 = E11/bla
        self.Q12 = nu12*E22/bla
        self.Q22 = E22/bla
        self.Q66 = G12

        # Ui: material invariants
        self.U1 = (1/8)*(3*self.Q11 + 3*self.Q22 + 2*self.Q12 + 4*self.Q66)
        self.U2 = (1/2)*(self.Q11 - self.Q22)
        self.U3 = (1/8)*(self.Q11 + self.Q22 - 2*self.Q12 - 4*self.Q66)
        self.U4 = (1/8)*(self.Q11 + self.Q22 + 6*self.Q12 - 4*self.Q66)
        self.U5 = (1/8)*(self.Q11 + self.Q22 - 2*self.Q12 + 4*self.Q66)

    def __repr__(self):
        " Display material properties"
        return f"""
Elastic modulus in the fibre direction : {self.E11:.2E} (Pa)
Elastic modulus in the transverse direction : {self.E22:.2E} (Pa)
In-plane shear modulus : {self.G12:.2E} (Pa)
Poisson's ratio relating transverse deformation and axial loading : {self.nu12}
Poisson's ratio relating axial deformation and transverse density : {self.nu21}
Volumic density : {self.density_volume} (g/m3)
Surface density : {self.density_area} (g/m2)
Ply thcikness : {self.ply_t} (m)

Material invariants:
U1 : {self.U1:.2E}
U2 : {self.U2:.2E}
U3 : {self.U3:.2E}
U4 : {self.U4:.2E}
U5 : {self.U5:.2E}

Ply stiffnesses:
Q11 : {self.Q11:.2E}
Q22 : {self.Q22:.2E}
Q12 : {self.Q12:.2E}
Q66 : {self.Q66:.2E}
"""

class MaterialDefinitionError(Exception):
    " Exceptions during the instanciation of a material"

if __name__ == "__main__":
    mat = Material()
    print(mat)
