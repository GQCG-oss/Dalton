module io_channels
  integer, parameter :: luinp   = 5
  integer, parameter :: lupri   = 6
  integer, parameter :: luerr   = 0
end module io_channels

module periodic_table
!
! Elemnt names are taken from the IUPAC recommendations, thus
! aluminium, caesium, and sulfur.
!
! Only the 112 "named" elements are included, no Uux-type names
! or symbols.
!

  integer, parameter       :: MAX_ELTS=112

  character(LEN=2)         :: elements(MAX_ELTS)
      data elements &
           /'H ',                                     'He', &
            'Li', 'Be' ,'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
            'Na', 'Mg' ,'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
            'K ', 'Ca' , &
            'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni' ,'Cu', 'Zn', &
                        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
            'Rb', 'Sr', &
            'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
                        'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
            'Cs', 'Ba', &
                'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', &
                'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
            'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
                        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
            'Fr', 'Ra', &
                 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', &
                 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', &
             'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'/
                

  integer, private        :: element_counter
  character(LEN=13)       :: element_names(MAX_ELTS)
  real (kind(0.d0))       :: atomic_weights(MAX_ELTS)

  data (element_names(element_counter), element_counter=1,10) &
      /'Hydrogen', 'Helium',   'Lithium', 'Beryllium', 'Boron', &
       'Carbon',   'Nitrogen', 'Oxygen',  'Fluorine',  'Neon'/
  data (atomic_weights(element_counter), element_counter=1,10) &
     / 1.0078246d0, 4.002601d0, 7.01600d0, 9.01218d0, 11.009307d0, &
       12.000000d0, 14.0030738d0, 15.9949141d0, 18.9984022d0, 19.992441d0/

  data (element_names(element_counter), element_counter=11,20) &
      /'Sodium', 'Magnesium', 'Aluminium', 'Silicon',   'Phosphorus', &
       'Sulfur', 'Chlorine',  'Argon',     'Potassium', 'Calcium'/
  data (atomic_weights(element_counter), element_counter=11,20) &
     / 22.9898d0, 23.98504d0, 26.98153d0, 27.976929d0, 30.973764d0, &
       31.9720727d0, 34.9688531d0, 39.962386d0, 38.96371d0, 39.96259d0/

  data (element_names(element_counter), element_counter=21,30) &
     /'Scandium', 'Titanium', 'Vanadium', 'Chromium', 'Manganese', &
      'Iron',     'Cobalt',   'Nickel',   'Copper',   'Zinc'/
  data (atomic_weights(element_counter), element_counter=21,30) &
     / 44.95592d0, 48.d0, 50.9440d0, 51.9405d0, 54.9380d0, &
       55.9349d0, 58.9332d0, 57.9353d0, 62.9296d0, 63.9291d0/

  data (element_names(element_counter), element_counter=31,40) &
     /'Gallium', 'Germanium', 'Arsenic',   'Selenium', 'Bromine', &
      'Krypton', 'Rubidium',  'Strontium', 'Yttrium',  'Zirconium'/
  data (atomic_weights(element_counter), element_counter=31,40) &
     / 68.9257d0, 73.9219d0, 74.9216d0, 79.9165d0, 78.91839d0, &
       83.91151d0, 84.9117d0, 87.9056d0, 88.9059d0, 89.9043d0/

  data (element_names(element_counter), element_counter=41,50) &
     /'Niobium',   'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium', &
      'Palladium', 'Silver',     'Cadmium',    'Indium',    'Tin'/
  data (atomic_weights(element_counter), element_counter=41,50) &
     / 92.9060d0, 97.9055d0, 98.d0, 101.9037d0, 102.9048d0, &
       107.90389d0, 106.90509d0, 113.9036d0, 114.9041d0, 120.d0/

  data (element_names(element_counter), element_counter=51,60) &
     /'Antimony', 'Tellurium', 'Iodine', 'Xenon',        'Caesium', &
      'Barium',   'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium'/
  data (atomic_weights(element_counter), element_counter=51,60) &
     / 120.9038d0, 129.9067d0, 126.90466d0, 131.90416d0, 132.9051d0, &
       137.9050d0, 138.9061d0, 139.9053d0, 140.9074d0, 141.9075d0/

  data (element_names(element_counter), element_counter=61,70) &
     /'Promethium', 'Samarium', 'Europium', 'Gadolinium', 'Terbium', &
      'Dysprosium', 'Holmium',  'Erbium',   'Thulium',    'Ytterbium'/
  data (atomic_weights(element_counter), element_counter=61,70) &
     / 145.d0, 151.9195d0, 152.9209d0, 157.9241d0, 159.9250d0, &
       163.9288d0, 164.9303d0, 165.9304d0, 168.9344d0, 173.9390d0/

  data (element_names(element_counter), element_counter=71,80) &
     /'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium', &
      'Osmium',   'Iridium', 'Platinum', 'Gold',     'Mercury'/
  data (atomic_weights(element_counter), element_counter=71,80) &
     / 174.9409d0, 179.9468d0, 180.9480d0, 183.9510d0, 186.9560d0, &
       192.d0, 192.9633d0, 194.9648d0, 196.9666d0, 201.970625d0/

  data (element_names(element_counter), element_counter=81,90) &
     /'Thallium', 'Lead',     'Bismuth', 'Polonium', 'Astatine', &
      'Radon',    'Francium', 'Radium',  'Actinium', 'Thorium'/
  data (atomic_weights(element_counter), element_counter=81,90) &
     / 204.9745d0, 207.9766d0, 208.9804d0, 209.d0, 210.d0, &
       222.d0, 223.d0, 226.d0, 227.d0, 232.d0/

  data (element_names(element_counter), element_counter=91,100) &
     /'Protactinium', 'Uranium',   'Neptunium',   'Plutonium',   'Americium', &
      'Curium',       'Berkelium', 'Californium', 'Einsteinium', 'Fermium'/
  data (atomic_weights(element_counter), element_counter=91,100) &
     / 231.d0, 238.d0, 237.d0, 244.d0, 243.d0, &
       247.d0, 247.d0, 251.d0, 252.d0, 257.d0/

  data (element_names(element_counter), element_counter=101,110) &
     /'Mendelevium', 'Nobelium', 'Lawrencium', 'Rutherfordium', 'Dubnium', &
      'Seaborgium',  'Bohrium',  'Hassium',    'Meitnerium',    'Darmstadtium'/
  data (atomic_weights(element_counter), element_counter=101,110) &
     / 258.d0, 259.d0, 262.d0, 267.d0, 268.d0, &
       269.d0, 270.d0, 269.d0, 278.d0, 281.d0/
       
  data (element_names(element_counter), element_counter=111,112) &
     /'Roentgenium', 'Copernicium'/
  data (atomic_weights(element_counter), element_counter=111,112) &
    / 281.d0, 285.d0/

end module periodic_table

module angular_momenta
  integer, parameter     ::  MAX_ANG_MOM=7

  character(LEN=1)       ::  l_labels(MAX_ANG_MOM)
  data  l_labels /'S', 'P', 'D', 'F', 'G', 'H', 'I' /

end module angular_momenta
