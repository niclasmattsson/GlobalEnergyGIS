export europe8, eurasia38, scand3, scand4, eurasia21, china6, europe13, europe21, europe22, europe54, syntheticdemandregions,
        testreg, caspian, NUTS_Europe, eurasia92

const scand3 = [
    "SWE"   GADM("Sweden")
    "NOR"   GADM("Norway")
    "FIN"   GADM("Finland","Åland")
]

const scand4 = [
    scand3;
    "DEN"   GADM("Denmark")
]

# const spainsubregions = subregions(GADM, "Spain")[[1:13; 15:18]]  # all subregions except Canary Islands
# const portugalsubregions = subregions(GADM, "Portugal")[[1; 3:11; 13:20]]  # all subregions except Azores and Madeira
# const europe8 = [
#     # Europe: "NOR","FRA","GER","UK","MED","BAL","SPA","CEN"
#     "NOR"   GADM("Sweden","Norway","Denmark","Finland","Åland","Faroe Islands")
#     "FRA"   GADM("France","Monaco")
#     "GER"   GADM("Germany","Netherlands","Belgium","Luxembourg")
#     "UK"    GADM("United Kingdom","Ireland","Guernsey","Isle of Man","Jersey")
#     "MED"   GADM("Greece","Bulgaria","Italy","San Marino","Vatican City","Slovenia","Croatia","Bosnia and Herzegovina","Serbia","Montenegro","Kosovo","Albania","Macedonia","Malta")
#     "BAL"   GADM("Poland","Estonia","Latvia","Lithuania")
#     "SPA"   (GADM(["Spain"], spainsubregions...), GADM(["Portugal"], portugalsubregions...), GADM("Andorra","Gibraltar"))
#     "CEN"   GADM("Austria","Switzerland","Czech Republic","Hungary","Slovakia","Romania","Liechtenstein")
# ]

const eurasia38 = [
    # Europe: "NOR","FRA","GER","UK","MED","BAL","SPA","CEN"
    "NOR"   GADM("Sweden","Norway","Denmark","Finland","Åland","Faroe Islands")
    "FRA"   GADM("France","Monaco")
    "GER"   GADM("Germany","Netherlands","Belgium","Luxembourg")
    "UK"    GADM("United Kingdom","Ireland","Guernsey","Isle of Man","Jersey")
    "MED"   GADM("Greece","Bulgaria","Italy","San Marino","Vatican City","Slovenia","Croatia","Bosnia and Herzegovina","Serbia","Montenegro","Kosovo","Albania","Macedonia","Malta")
    "BAL"   GADM("Poland","Estonia","Latvia","Lithuania")
    "SPA"   GADM("Spain","Portugal","Andorra","Gibraltar")
    "CEN"   GADM("Austria","Switzerland","Czech Republic","Hungary","Slovakia","Romania","Liechtenstein")

    # Eastern Europe & Middle East: "Belarus/Ukraine/Moldova","Turkey/Caucasus","Middle East","Iran","Arabian peninsula"
    "BUK"   GADM("Belarus","Ukraine","Moldova")
    "TCC"   GADM("Turkey","Georgia","Armenia","Azerbaijan")
    "MEA"   GADM("Israel","Palestina","Lebanon","Jordan","Syria","Iraq","Kuwait")
    "IRN"   GADM("Iran")
    "ARB"   GADM("Saudi Arabia","Yemen","Oman","United Arab Emirates","Bahrain","Qatar")
    
    # Central & South-Central Asia: "Central Asia","South-Central Asia"
    "KZK"   GADM("Kazakhstan")
    "CAS"   GADM("Uzbekistan","Turkmenistan","Tajikistan","Kyrgyzstan")
    "SCA"   GADM("Afghanistan","Pakistan")
    
    # India: "North","West","Central","South","East","Northeast"
    # https://en.wikipedia.org/wiki/Administrative_divisions_of_India
    # Excluding island groups: "Andaman and Nicobar Islands","Lakshadweep"
    # South also includes Sri Lanka
    # Northeast also includes Nepal,Bhutan and Bangladesh
    "IN_N"  GADM(["India"], "Jammu and Kashmir","Himachal Pradesh","Punjab","Rajasthan","Chandigarh","Haryana","NCT of Delhi","Uttarakhand","Uttar Pradesh")
    "IN_W"  GADM(["India"], "Gujarat","Goa","Maharashtra","Dadra and Nagar Haveli","Daman and Diu")
    "IN_C"  GADM(["India"], "Madhya Pradesh","Chhattisgarh")
    "IN_S"  (GADM(["India"], "Karnataka","Kerala","Tamil Nadu","Andhra Pradesh","Telangana","Puducherry"), GADM("Sri Lanka"))
    "IN_E"  GADM(["India"], "Odisha","Jharkhand","West Bengal","Bihar")
    "IN_NE" (GADM(["India"], "Sikkim","Assam","Meghalaya","Tripura","Mizoram","Manipur","Nagaland","Arunachal Pradesh"),
                GADM("Nepal","Bhutan","Bangladesh"))

    # Southeast Asia: "Central Asia","South-Central Asia"
    # also includes Peninsular Malaysia,https://en.wikipedia.org/wiki/Peninsular_Malaysia
    "SEA"   (GADM("Myanmar","Thailand","Laos","Vietnam","Cambodia","Singapore"),
                GADM(["Malaysia"], "Perlis","Kedah","Pulau Pinang","Perak","Kelantan","Trengganu","Pahang","Selangor","Kuala Lumpur","Putrajaya","Negeri Sembilan","Melaka","Johor"))
    
    # Russia: "Northwest","Central","Southwest","Volga","Ural","Siberia","East"
    # https://en.wikipedia.org/wiki/Federal_districts_of_Russia
    # questionable: Novgorod,Altay,Yevrey,Maga Buryatdan/Magadan
    "RU_NW"  GADM(["Russia"], "Arkhangel'sk","Vologda","Kaliningrad","Karelia","Komi","Leningrad","Murmansk","Nenets","Novgorod","Pskov","City of St. Petersburg")
    "RU_C"   GADM(["Russia"], "Belgorod","Bryansk","Vladimir","Voronezh","Ivanovo","Kaluga","Kostroma","Kursk","Lipetsk","Moscow City","Moskva","Orel","Ryazan'","Smolensk","Tambov","Tver'","Tula","Yaroslavl'")
    "RU_SW"  GADM(["Russia"], "Adygey","Astrakhan'","Volgograd","Kalmyk","Krasnodar","Rostov","Dagestan","Ingush","Kabardin-Balkar","Karachay-Cherkess","North Ossetia","Stavropol'","Chechnya")
    "RU_VL"  GADM(["Russia"], "Bashkortostan","Kirov","Mariy-El","Mordovia","Nizhegorod","Orenburg","Penza","Perm'","Samara","Saratov","Tatarstan","Udmurt","Ul'yanovsk","Chuvash")
    "RU_UR"  GADM(["Russia"], "Kurgan","Sverdlovsk","Tyumen'","Khanty-Mansiy","Chelyabinsk","Yamal-Nenets")
    "RU_SB"  GADM(["Russia"], "Altay","Gorno-Altay","Irkutsk","Kemerovo","Krasnoyarsk","Novosibirsk","Omsk","Tomsk","Tuva","Khakass")
    "RU_E"   GADM(["Russia"], "Amur","Buryat","Yevrey","Zabaykal'ye","Kamchatka","Maga Buryatdan","Primor'ye","Sakha","Sakhalin","Khabarovsk","Chukot")
    
    # China: "North" (Huáběi),"Northeast" (Dōngběi),"East" (Huádōng),"South Central" (Zhōngnán),"Southwest" (Xīnán),"Northwest" (Xīběi)
    # https://en.wikipedia.org/wiki/List_of_regions_of_China
    # East also includes Taiwan
    # South Central also includes Hong Kong and Macao
    "CH_N"   GADM(["China"], "Beijing","Tianjin","Hebei","Shanxi","Nei Mongol")
    "CH_NE"  GADM(["China"], "Liaoning","Jilin","Heilongjiang")
    "CH_E"   (GADM(["China"], "Shanghai","Jiangsu","Zhejiang","Anhui","Fujian","Jiangxi","Shandong"), GADM("Taiwan"))
    "CH_SC"  (GADM(["China"], "Henan","Hubei","Hunan","Guangdong","Guangxi","Hainan"), GADM("Hong Kong","Macao"))
    "CH_SW"  GADM(["China"], "Chongqing","Sichuan","Guizhou","Yunnan","Xizang")
    "CH_NW"  GADM(["China"], "Shaanxi","Gansu","Qinghai","Ningxia Hui","Xinjiang Uygur")

    # Other
    "MON"    GADM("Mongolia")
    "JKR"    GADM("Japan","South Korea","North Korea")
]

const eurasia21 = eurasia38[[1:10; 14:15; 25:27; 31:36], :]
const china6 = eurasia38[31:36, :]

const europe13 = [
    # BNL = Benelux, BTC = Baltic, BKN = Balkan, CEN = Central
    "SWE"   GADM("Sweden")
    "NOR"   GADM("Norway")
    "DEN"   GADM("Denmark","Faroe Islands")
    "FIN"   GADM("Finland","Åland")
    "BNL"   GADM("Netherlands","Belgium","Luxembourg")
    "GER"   GADM("Germany")
    "FRA"   GADM("France","Monaco")
    "UK"    GADM("United Kingdom","Ireland","Guernsey","Isle of Man","Jersey")
    "ITA"   GADM("Italy","San Marino","Vatican City")
    "BKN"   GADM("Greece","Bulgaria","Romania","Slovenia","Croatia","Bosnia and Herzegovina","Serbia","Montenegro","Kosovo","Albania","Macedonia","Malta","Moldova")
    "BTC"   GADM("Poland","Estonia","Latvia","Lithuania")
    "SPA"   GADM("Spain","Portugal","Andorra","Gibraltar")
    "CEN"   GADM("Austria","Switzerland","Czech Republic","Hungary","Slovakia","Liechtenstein")
]

# Wikipedia's definition of Europe
const europe21 = [
    europe13;
    "BLR"   GADM("Belarus")
    "UKR"   GADM("Ukraine")
    "TUR"   GADM("Turkey")
    "CCS"   GADM("Georgia","Armenia","Azerbaijan")
    eurasia38[24:27, :]
]

# Regions for the ELIN/EPOD models
const europe54 = [
    "AT"    NUTS("AT")
    "BE"    NUTS("BE")
    "BG"    NUTS("BG")
    "CY"    NUTS("CY")
    "CZ"    NUTS("CZ")
    "DE1"   NUTS("DE11","DE12","DE13","DE14","DE21","DE22","DE23","DE24","DE25","DE26","DE27")
    "DE2"   NUTS("DEB1","DEB2","DEB3","DEC0")
    "DE3"   NUTS("DE71","DE72","DE73","DE91","DE92","DEA1","DEA2","DEA3","DEA4","DEA5","DEG0")
    "DE4"   NUTS("DE30","DE40","DE60","DE80","DED4","DED2","DED5","DEF0","DEE0")
    "DE5"   NUTS("DE50","DE93","DE94")
    "DK1"   NUTS("DK01","DK02")
    "DK2"   NUTS("DK03","DK04","DK05")
    "EE"    NUTS("EE")
    "FI"    NUTS("FI")                # Åland is included in NUTS("FI")  
    "GR"    NUTS("EL")                # ΕΛΛΑΔΑ (Greece in Greek)   
    "HU"    NUTS("HU")
    "IE"    NUTS("IE")
    "LV"    NUTS("LV")
    "LT"    NUTS("LT")
    "LU"    NUTS("LU")
    "MT"    NUTS("MT")
    "NL"    NUTS("NL")
    "PT"    NUTS("PT11","PT15","PT16","PT17","PT18")    # excludes small islands
    "RO"    NUTS("RO")
    "SI"    NUTS("SI")
    "SK"    NUTS("SK")
    "CH"    NUTS("CH","LI")                   # includes Liechtenstein
    "IS"    NUTS("IS")
 #   "BO"    GADM("Bosnia and Herzegovina")    # Bosnia and Herzegovina completely missing in NUTS
    "CR"    NUTS("HR")                        # Hrvatska (Croatia in Croatian)
 #   "MC"    NUTS("MK")                        # Macedonia
    "ES1"   NUTS("ES11","ES12","ES13")  
    "ES2"   NUTS("ES21","ES22","ES23","ES24","ES51","ES52","ES53","ES62")
    "ES3"   NUTS("ES30","ES41","ES42")
    "ES4"   NUTS("ES61","ES43")
    "FR1"   NUTS("FRJ1","FRL0","FRM0")
    "FR2"   NUTS("FRI1","FRJ2","FRI2")
    "FR3"   NUTS("FRG0","FRH0","FRI3")
    "FR4"   NUTS("FRF3","FRF1","FRC2")
    "FR5"   NUTS("FR10","FRF2","FRE2","FRD2","FRB0","FRD1","FRC1","FRE1","FRK2","FRK1")
    "IT1"   NUTS("ITC","ITH")
    "IT2"   NUTS("ITI","ITG2")
    "IT3"   NUTS("ITF","ITG1")
    "NO1"   NUTS("NO01","NO02","NO03","NO04","NO05")
    "NO2"   NUTS("NO06")
    "NO3"   NUTS("NO07")
    "PO1"   NUTS("PL21","PL22","PL51","PL52")
    "PO2"   NUTS("PL71","PL72","PL81","PL82","PL84","PL91","PL92")
    "PO3"   NUTS("PL41","PL42","PL43","PL61","PL62","PL63")
    "SE1"   NUTS("SE22")
    "SE2"   NUTS("SE11","SE12","SE21","SE23","SE31")
    "SE3"   NUTS("SE32")
    "SE4"   NUTS("SE33")
    "UK1"   NUTS("UKC1","UKC2","UKD1","UKD3","UKD4","UKD5","UKD6","UKD7","UKL1","UKL2","UKG1","UKG2","UKG3",
                 "UKE1","UKE2","UKE3","UKE4","UKF1","UKF2","UKF3","UKH1","UKH2","UKH3","UKI3","UKI4","UKI5","UKI6","UKI7",
                 "UKK1","UKK2","UKK3","UKK4","UKJ1","UKJ2","UKJ3","UKJ4")
    "UK2"   NUTS("UKM5","UKM6","UKM7","UKM8","UKM9")
    "UK3"   NUTS("UKN0", "UKN1")
]

# https://ec.europa.eu/eurostat/statistics-explained/index.php/Glossary:EU_enlargements
# https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Candidate_countries
const NUTS_Europe = [
    "EU12"  NUTS("BE","DK","FR","DE","EL","IE","IT","LU","NL","PT","ES","UK")
    "EU15"  NUTS("AT","FI","SE")
    "EFTA"  NUTS("IS","LI","NO","CH")
    "EU25"  NUTS("CY","CZ","EE","HU","LV","LT","MT","PL","SK","SI")
    "EU28"  NUTS("BG","RO","HR")
    "CAND"  NUTS("ME","MK","AL","RS","TR","BA","XK")
]

# Non NUTS regions and neighboring countries to European regions (onshore or offshore)
# https://ec.europa.eu/eurostat/statistical-atlas/gis/viewer/?config=typologies.json
const non_NUTS_Europe = [
    "mainland"  GADM("Monaco", "Andorra", "San Marino", "Vatican City")
    "islands"   GADM("Faroe Islands", "Isle of Man", "Guernsey", "Jersey", "Svalbard and Jan Mayen")
    "east"      GADM("Moldova", "Ukraine", "Belarus", "Russia")
    "south"     GADM("Morocco", "Algeria", "Tunisia", "Libya", "Egypt", "Israel", "Jordan", "Lebanon", "Syria")
    "west"      GADM("Greenland")
]

const europe22 = [
    "SE1"   NUTS("SE332")
    "SE2"   NUTS("SE313", "SE321", "SE322", "SE331")
    "SE3"   NUTS("SE110", "SE121", "SE122", "SE123", "SE124", "SE125", "SE211", "SE213", "SE214", "SE232", "SE311", "SE312")
    "SE4"   NUTS("SE212", "SE221", "SE224", "SE231")
    "NO1"   NUTS("NO011", "NO012", "NO021", "NO022", "NO031")
    "NO2"   NUTS("NO033", "NO034", "NO041", "NO042", "NO043")
    "NO3"   NUTS("NO053", "NO060")
    "NO4"   NUTS("NO071", "NO072", "NO073")
    "NO5"   NUTS("NO032", "NO051", "NO052")
    "DEN"   NUTS("DK")
    "FI1"   NUTS("FI193", "FI194", "FI195", "FI196", "FI197", "FI1B1", "FI1C1", "FI1C2", "FI1C3", "FI1C4", "FI1C5",
                    "FI1D1", "FI1D2", "FI1D3", "FI1D5", "FI200")
    "FI2"   NUTS("FI1D7", "FI1D8", "FI1D9")
    "BNL"   NUTS("NL","BE","LU")
    "GER"   NUTS("DE")
    "POL"   NUTS("PL")
    "BAL"   NUTS("EE", "LV", "LT")
    "IRL"   NUTS("IE")
    "UK"    NUTS("UK")
    "FRA"   (NUTS("FR"), GADM("Monaco"))
    "CEN"   NUTS("AT", "CH", "LI", "CZ", "HU", "SK", "RO")  # Austria, Switzerland, Liechtenstein, Czech Rep, Hungary, Slovakia, Romania
    "SPA"   (NUTS("ES", "PT"), GADM("Andorra"))             # Spain includes Gibraltar
    "MED"   (NUTS("IT", "SI", "HR", "RS", "BG", "ME",       # Italy, Slovenia, Croatia, Serbia, Bularia, Montenegro
                  "AL", "MK", "EL"),                        # Albania, North Macedonia, Greece
            GADM("San Marino", "Vatican City", "Bosnia and Herzegovina", "Kosovo"))       # these are not included in NUTS-2016
    # Not included in any region:  Iceland, Faroe Islands, Cyprus, Malta, Isle of Man, Guernsey, Jersey, Jan Mayen, Svalbard, Bear Island,
    #                               Moldova, Ukraine, Belarus, Kaliningrad (or any part of Russia)
    # Islands that ARE included:  Canarias, Mallorca and other Balearic islands (ES); Madeira, Azores (PT); Åland (FI);
    #                               Corsica, Reunion, Mayotte, Guadeloupe, Martinique, Guyana (FR); Sicily, Sardinia (IT);
    #                               Orkney Islands, Shetney Islands, Western Isles, Isle of Wight (UK); Crete (EL)
    # But note that non-European territories are eliminated when the region bounding box is set to [34 -11; 72 32].
]

const france13 = [
    "Auvergne-Rhône-Alpes"         GADM(["France"], "Auvergne-Rhône-Alpes")
    "Bourgogne-Franche-Comté"      GADM(["France"], "Bourgogne-Franche-Comté")
    "Bretagne"                     GADM(["France"], "Bretagne")
    "Centre-Val de Loire"          GADM(["France"], "Centre-Val de Loire")
    "Corse"                        GADM(["France"], "Corse")
    "Grand Est"                    GADM(["France"], "Grand Est")
    "Hauts-de-France"              GADM(["France"], "Hauts-de-France")
    "Île-de-France"                GADM(["France"], "Île-de-France")
    "Normandie"                    GADM(["France"], "Normandie")
    "Nouvelle-Aquitaine"           GADM(["France"], "Nouvelle-Aquitaine")
    "Occitanie"                    GADM(["France"], "Occitanie")
    "Pays de la Loire"             GADM(["France"], "Pays de la Loire")
    "Provence-Alpes-Côte d'Azur"   GADM(["France"], "Provence-Alpes-Côte d'Azur")
]

const usa9 = [
    "New England"       GADM(["United States"], "Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont")
    "Mid-Atlantic"      GADM(["United States"], "New Jersey", "New York", "Pennsylvania")
    "Central NE"        GADM(["United States"], "Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin")
    "Central NW"        GADM(["United States"], "Iowa", "Kansas", "Minnesota", "Missouri", "Nebraska", "North Dakota", "South Dakota")
    "South Atlantic"    GADM(["United States"], "Delaware", "Florida", "Georgia", "Maryland", "North Carolina", "South Carolina", "Virginia", "District of Columbia", "West Virginia")
    "Central SE"        GADM(["United States"], "Alabama", "Kentucky", "Mississippi", "Tennessee")
    "Central SW"        GADM(["United States"], "Arkansas", "Louisiana", "Oklahoma", "Texas")
    "Mountain"          GADM(["United States"], "Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico", "Utah", "Wyoming")
    "Pacific"           GADM(["United States"], "California", "Oregon", "Washington") # "Alaska", "Hawaii",
]

const denmark5 = [
    "Hovedstaden"       GADM(["Denmark"], "Hovedstaden")
    "Midtjylland"       GADM(["Denmark"], "Midtjylland")
    "Nordjylland"       GADM(["Denmark"], "Nordjylland")
    "Sjælland"          GADM(["Denmark"], "Sjælland")
    "Syddanmark"        GADM(["Denmark"], "Syddanmark")
]

# Test region to see if mixed GADM/NUTS notation works
const testreg = [
    "NO"     GADM("Norway")
    "DK"     NUTS("DK")
    "ISL"    (NUTS("SE214"), GADM("Åland"), GADM(["Sweden","Kalmar"], "Borgholm","Mörbylånga"))  # Gotland, Åland, Öland
]

const caspian = [
    "RU_SW"  GADM(["Russia"], "Astrakhan'","Kalmyk", "Dagestan")
    "KZK"    GADM(["Kazakhstan"], "Atyrau", "Mangghystau")
    "AZR"    GADM("Azerbaijan")
    "IRN"    GADM(["Iran"], "Gilan", "Golestan", "Mazandaran")
    "TRK"    GADM(["Turkmenistan"], "Balkan")
]

const mapregions = [
    "FRA"   GADM("France","Monaco")
    "KZK"   GADM("Kazakhstan")
    "CH_E"   (GADM(["China"], "Shanghai","Jiangsu","Zhejiang","Anhui","Fujian","Jiangxi","Shandong"), GADM("Taiwan"))    
]

const syntheticdemandregions = [
    "Argentina", "Australia", "Austria", "Belgium", "Bosnia and Herzegovina", "Brazil", "Bulgaria", "Canada", "Chile",
    "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary",
    "Iceland", "Ireland", "Italy", "Kenya", "Latvia", "Lithuania", "Macedonia", "Mexico", "Netherlands", "New Zealand",
    "Norway", "Poland", "Portugal", "Romania", "Saudi Arabia", "Serbia", "Slovakia", "Slovenia", "South Korea", "Spain",
    "Sri Lanka", "Sweden", "Switzerland", "Turkey", "United Kingdom"
]

const hovedstaden15 = [
    "København C"     GADM(["Denmark", "Hovedstaden"], "København","Frederiksberg","Tårnby","Dragør")
    "København W"     GADM(["Denmark", "Hovedstaden"], "Albertslund","Glostrup","Rødovre","Vallensbæk","Brøndby","Hvidovre")
    "København N"     GADM(["Denmark", "Hovedstaden"], "Lyngby-Taarbæk","Gladsaxe","Gentofte")
    "København NW"    GADM(["Denmark", "Hovedstaden"], "Furesø","Ballerup","Herlev")
    "Hovedstaden SW"  GADM(["Denmark", "Hovedstaden"], "Høje Taastrup", "Ishøj")
    "Allerød"         GADM(["Denmark", "Hovedstaden"], "Allerød")
    "Bornholm"        GADM(["Denmark", "Hovedstaden"], "Bornholm")  # Christiansø
    "Egedal"          GADM(["Denmark", "Hovedstaden"], "Egedal")
    "Helsingør"       GADM(["Denmark", "Hovedstaden"], "Fredensborg","Helsingør")
    "Frederikssund"   GADM(["Denmark", "Hovedstaden"], "Frederikssund")
    "Gribskov"        GADM(["Denmark", "Hovedstaden"], "Gribskov")
    "Halsnæs"         GADM(["Denmark", "Hovedstaden"], "Halsnæs")
    "Hillerød"        GADM(["Denmark", "Hovedstaden"], "Hillerød")
    "Hørsholm"        GADM(["Denmark", "Hovedstaden"], "Hørsholm")
    "Rudersdal"       GADM(["Denmark", "Hovedstaden"], "Rudersdal")
]

# no small islands included
const eurasia92 = [
    # Europe: "NOR","FRA","GER","UK","MED","BAL","SPA","CEN"
    "NOR_SEN"   GADM(["Sweden"], "Norrbotten","Västerbotten","Västernorrland","Jämtland","Gävleborg","Dalarna","Värmland")
    "NOR_SES"   GADM(["Sweden"], "Uppsala","Västmanland","Orebro","Stockholm","Södermanland","Västra Götaland","Östergötland","Jönköping","Kalmar","Gotland","Halland","Kronoberg","Skåne","Blekinge")
    "NOR_NO"    GADM("Norway")
    "NOR_DK"    GADM("Denmark")
    "NOR_FI"    GADM("Finland","Åland")
    "FRA_N"   GADM(["France"], "Hauts-de-France","Normandie")
    "FRA_W"   GADM(["France"], "Bretagne","Pays de la Loire","Nouvelle-Aquitaine")
    "FRA_C"   GADM(["France"], "Île-de-France","Grand Est","Centre-Val de Loire","Bourgogne-Franche-Comté")
    "FRA_S"   (GADM(["France"], "Auvergne-Rhône-Alpes","Occitanie","Provence-Alpes-Côte d'Azur"), GADM("Monaco"))
    "GER_BNL" GADM("Netherlands","Belgium","Luxembourg")
    "GER_N"   GADM(["Germany"], "Schleswig-Holstein","Hamburg","Mecklenburg-Vorpommern","Bremen","Niedersachsen")
    "GER_W"   GADM(["Germany"], "Nordrhein-Westfalen","Hessen","Rheinland-Pfalz","Saarland")
    "GER_E"   GADM(["Germany"], "Brandenburg","Berlin","Sachsen-Anhalt","Thüringen","Sachsen")
    "GER_S"   GADM(["Germany"], "Baden-Württemberg","Bayern")
    "UK_ESW"  GADM(["United Kingdom"], "England","Scotland","Wales")
    "UK_I"    (GADM("Ireland"), GADM(["United Kingdom"], "Northern Ireland"))
    "MED_IT"   GADM("Italy","San Marino","Vatican City")
    "MED_BK"   GADM("Slovenia","Croatia","Bosnia and Herzegovina","Serbia","Montenegro","Kosovo")
    "MED_GR"   GADM("Greece","Albania","Macedonia")
    "BAL_PL"   GADM("Poland")
    "BAL_ELL"  (GADM("Estonia","Latvia","Lithuania"), GADM(["Russia"], "Kaliningrad"))
    "SPA_NW"   GADM(["Spain"], "Galicia","Principado de Asturias","Cantabria","País Vasco","Comunidad Foral de Navarra","La Rioja","Castilla y León")
    "SPA_C"    GADM(["Spain"], "Aragón","Castilla-La Mancha","Comunidad de Madrid","Extremadura")
    "SPA_SE"   (GADM(["Spain"], "Cataluña", "Comunidad Valenciana", "Región de Murcia", "Andalucía"), GADM("Andorra","Gibraltar"))
    "SPA_PT"   GADM("Portugal")
    "CEN_AT"   GADM("Austria")
    "CEN_CH"   GADM("Switzerland","Liechtenstein")
    "CEN_CZ"   GADM("Czech Republic")
    "CEN_HU"   GADM("Hungary")
    "CEN_SK"   GADM("Slovakia")
    "CEN_RO"   GADM("Romania","Moldova")
    "CEN_BG"   GADM("Bulgaria")

    # Eastern Europe & Middle East: "Belarus/Ukraine/Moldova","Turkey/Caucasus","Middle East","Iran","Arabian peninsula"
    "BUK_BY"   GADM("Belarus")
    "BUK_UA"   GADM("Ukraine")
    "TCC_TR"   GADM("Turkey")
    "TCC_CC"   GADM("Georgia","Armenia","Azerbaijan")
    "MEA_SY"   GADM("Syria")
    "MEA_IPL"  GADM("Israel","Palestina","Lebanon")
    "MEA_JO"   GADM("Jordan")
    "MEA_IQ"   GADM("Iraq","Kuwait")
    "ARB_SAN"  GADM(["Saudi Arabia"], "Tabuk","Al Jawf","Al Madinah","Ha'il","Al Quassim","Al Hudud ash Shamaliyah")
    "ARB_SAS"  GADM(["Saudi Arabia"], "Makkah","Ar Riyad","Al Bahah","`Asir","Jizan","Najran","Ash Sharqiyah")
    "ARB_S"   GADM("Yemen","Oman")
    "ARB_E"   GADM("United Arab Emirates","Bahrain","Qatar")
    "IRN_NW"  GADM(["Iran"], "Ardebil","East Azarbaijan","Gilan","Kordestan","West Azarbaijan","Zanjan")
    "IRN_W"   GADM(["Iran"], "Hamadan","Ilam","Kermanshah","Khuzestan","Lorestan","Markazi")
    "IRN_N"   GADM(["Iran"], "Alborz","Golestan","Mazandaran","Qazvin","Qom","Semnan","Tehran")
    "IRN_S"   GADM(["Iran"], "Bushehr","Chahar Mahall and Bakhtiari","Fars","Hormozgan","Esfahan","Kohgiluyeh and Buyer Ahmad")
    "IRN_E"   GADM(["Iran"], "Kerman","North Khorasan","Razavi Khorasan","Sistan and Baluchestan","South Khorasan","Yazd")

    # Central & South-Central Asia: "Central Asia","South-Central Asia"
    "KZK_E"   GADM(["Kazakhstan"], "Almaty","East Kazakhstan")
    "KZK_N"   GADM(["Kazakhstan"], "Qostanay","North Kazakhstan","Aqmola","Pavlodar")
    "KZK_W"   GADM(["Kazakhstan"], "West Kazakhstan","Atyrau","Aqtöbe","Mangghystau")
    "KZK_SC"  GADM(["Kazakhstan"], "Qaraghandy","Qyzylorda","South Kazakhstan","Zhambyl")
    "CAS_UZ"   GADM("Uzbekistan")
    "CAS_TM"   GADM("Turkmenistan")
    "CAS_KT"   GADM("Tajikistan","Kyrgyzstan")
    "SCA_AF"   GADM("Afghanistan")
    "SCA_PK"   GADM("Pakistan")
    
    # India: "North","West","Central","South","East","Northeast"
    # https://en.wikipedia.org/wiki/Administrative_divisions_of_India
    # Excluding island groups: "Andaman and Nicobar Islands","Lakshadweep"
    # South also includes Sri Lanka
    # Northeast also includes Nepal,Bhutan and Bangladesh
    "IN_N"  GADM(["India"], "Jammu and Kashmir","Himachal Pradesh","Punjab","Rajasthan","Chandigarh","Haryana","NCT of Delhi","Uttarakhand","Uttar Pradesh")
    "IN_W"  GADM(["India"], "Gujarat","Goa","Maharashtra","Dadra and Nagar Haveli","Daman and Diu")
    "IN_C"  GADM(["India"], "Madhya Pradesh","Chhattisgarh")
    "IN_S"  (GADM(["India"], "Karnataka","Kerala","Tamil Nadu","Andhra Pradesh","Telangana","Puducherry"), GADM("Sri Lanka"))
    "IN_E"  GADM(["India"], "Odisha","Jharkhand","West Bengal","Bihar","Sikkim")
    "IN_NE" GADM(["India"], "Assam","Meghalaya","Tripura","Mizoram","Manipur","Nagaland","Arunachal Pradesh")

    # Southeast Asia: "Central Asia","South-Central Asia"
    # also includes Peninsular Malaysia,https://en.wikipedia.org/wiki/Peninsular_Malaysia
    "SEA_NB"   GADM("Nepal","Bhutan")
    "SEA_BD"   GADM("Bangladesh")
    "SEA_VLC"  GADM("Laos","Vietnam","Cambodia")
    "SEA_MM"   GADM("Myanmar")
    "SEA_TH"   GADM("Thailand")
    "SEA_MY"   (GADM("Singapore"), GADM(["Malaysia"], "Perlis","Kedah","Pulau Pinang","Perak","Kelantan","Trengganu","Pahang","Selangor","Kuala Lumpur","Putrajaya","Negeri Sembilan","Melaka","Johor"))
    
    # Russia: "Northwest","Central","Southwest","Volga","Ural","Siberia","East"
    # https://en.wikipedia.org/wiki/Federal_districts_of_Russia
    # questionable: Novgorod,Altay,Yevrey,Maga Buryatdan/Magadan
    "RU_NW"  GADM(["Russia"], "Arkhangel'sk","Vologda","Karelia","Komi","Leningrad","Murmansk","Nenets","Novgorod","Pskov","City of St. Petersburg")
    "RU_C"   GADM(["Russia"], "Belgorod","Bryansk","Vladimir","Voronezh","Ivanovo","Kaluga","Kostroma","Kursk","Lipetsk","Moscow City","Moskva","Orel","Ryazan'","Smolensk","Tambov","Tver'","Tula","Yaroslavl'")
    "RU_SW"  GADM(["Russia"], "Adygey","Astrakhan'","Volgograd","Kalmyk","Krasnodar","Rostov","Dagestan","Ingush","Kabardin-Balkar","Karachay-Cherkess","North Ossetia","Stavropol'","Chechnya")
    "RU_VL"  GADM(["Russia"], "Bashkortostan","Kirov","Mariy-El","Mordovia","Nizhegorod","Orenburg","Penza","Perm'","Samara","Saratov","Tatarstan","Udmurt","Ul'yanovsk","Chuvash")
    "RU_UR"  GADM(["Russia"], "Kurgan","Sverdlovsk","Tyumen'","Khanty-Mansiy","Chelyabinsk","Yamal-Nenets")
    "RU_SB"  GADM(["Russia"], "Irkutsk","Krasnoyarsk","Tuva")
    "RU_S"  GADM(["Russia"], "Altay","Gorno-Altay","Kemerovo","Novosibirsk","Omsk","Tomsk","Khakass")
    "RU_E"   GADM(["Russia"], "Amur","Buryat","Yevrey","Zabaykal'ye","Kamchatka","Maga Buryatdan","Primor'ye","Sakha","Sakhalin","Khabarovsk","Chukot")
    
    # China: "North" (Huáběi),"Northeast" (Dōngběi),"East" (Huádōng),"South Central" (Zhōngnán),"Southwest" (Xīnán),"Northwest" (Xīběi)
    # https://en.wikipedia.org/wiki/List_of_regions_of_China
    # East also includes Taiwan
    # South Central also includes Hong Kong and Macao
    "CH_N"    GADM(["China"], "Nei Mongol")
    "CH_NE"   GADM(["China"], "Liaoning","Jilin","Heilongjiang")
    "CH_CNE"  GADM(["China"], "Beijing","Tianjin","Hebei","Shanxi","Shandong")
    "CH_CSE"  (GADM(["China"], "Shanghai","Jiangsu","Zhejiang","Anhui","Fujian","Jiangxi"), GADM("Taiwan"))
    "CH_CS"   GADM(["China"], "Henan","Hubei","Hunan")
    "CH_CSW"  GADM(["China"], "Chongqing","Sichuan","Guizhou","Yunnan")
    "CH_CNW"  GADM(["China"], "Shaanxi","Gansu","Qinghai","Ningxia Hui")
    "CH_S"    (GADM(["China"], "Guangdong","Guangxi","Hainan"), GADM("Hong Kong","Macao"))
    "CH_TB"   GADM(["China"], "Xizang")
    "CH_NW"   GADM(["China"], "Xinjiang Uygur")

    # Other
    "MON"      GADM("Mongolia")
    "JKR_JP"    GADM("Japan")
    "JKR_NK"    GADM("North Korea")
    "JKR_SK"    GADM("South Korea")
]
