export europe8, eurasia38, scand3, eurasia21, china6, europe13, europe56

const scand3 = [
    "SWE"   GADM("Sweden")
    "NOR"   GADM("Norway")
    "FIN"   GADM("Finland")
]

const europe8 = [
    # Europe: "NOR","FRA","GER","UK","MED","BAL","SPA","CEN"
    "NOR"   GADM("Sweden","Norway","Denmark","Finland","Åland","Faroe Islands")
    "FRA"   GADM("France","Monaco")
    "GER"   GADM("Germany","Netherlands","Belgium","Luxembourg")
    "UK"    GADM("United Kingdom","Ireland","Guernsey","Isle of Man","Jersey")
    "MED"   GADM("Greece","Bulgaria","Italy","San Marino","Vatican City","Slovenia","Croatia","Bosnia and Herzegovina","Serbia","Montenegro","Kosovo","Albania","Macedonia","Malta")
    "BAL"   GADM("Poland","Estonia","Latvia","Lithuania")
    "SPA"   GADM("Spain","Portugal","Andorra","Gibraltar")
    "CEN"   GADM("Austria","Switzerland","Czech Republic","Hungary","Slovakia","Romania","Liechtenstein")
]

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
    "BKN"   GADM("Greece","Bulgaria","Romania","Slovenia","Croatia","Bosnia and Herzegovina","Serbia","Montenegro","Kosovo","Albania","Macedonia","Malta")
    "BTC"   GADM("Poland","Estonia","Latvia","Lithuania")
    "SPA"   GADM("Spain","Portugal","Andorra","Gibraltar")
    "CEN"   GADM("Austria","Switzerland","Czech Republic","Hungary","Slovakia","Liechtenstein")
]

# Regions for the ELIN/EPOD models
const europe56 = [
    "AT"    NUTS("AT")
    "BE"    NUTS("BE")
    "BG"    NUTS("BG")
    "CY"    NUTS("CY")
    "CZ"    NUTS("CZ")
    "DE1"   NUTS("DE11","DE12","DE13","DE14","DE21","DE22","DE23","DE24","DE25","DE26","DE27")
    "DE2"   NUTS("DEB1","DEB2","DEB3","DEC0")      
    "DE3"   NUTS("DE71","DE72","DE73","DE91","DE92","DEA1","DEA2","DEA3","DEA4","DEA5","DEG0")
    "DE4"   NUTS("DE30","DE41","DE42","DE60","DE80","DED1","DED2","DED3","DEE0","DEF0")
    "DE5"   NUTS("DE50","DE93","DE94")
    "DK1"   NUTS("DK01","DK02")
    "DK2"   NUTS("DK03","DK04","DK05")
    "EE"    NUTS("EE")   
    "FI"    NUTS("FI")    # Åland is included in NUTS("FI")  
    "GR"    NUTS("EL")    # ΕΛΛΑΔΑ (Greece in Greek)   
    "HU"    NUTS("HU")
    "IE"    NUTS("IE")
    "LV"    NUTS("LV")
    "LT"    NUTS("LT")
    "LU"    NUTS("LU")
    "MT"    NUTS("MT")
    "NL"    NUTS("NL")   
    "PT"    NUTS("PT")
    "RO"    NUTS("RO")
    "SI"    NUTS("SI")
    "SK"    NUTS("SK")
    "CH"    NUTS("CH","LI")
    "IS"    NUTS("IS")    
    "BO"    NUTS("BO")    # Bosnia and Herzegovina completely missing in NUTS (?)
    "CR"    NUTS("HR")    # Hrvatska (Croatia in Croatian)
    "MC"    NUTS("MK")    # Macedonia
    "ES1"   NUTS("ES11","ES12","ES13")  
    "ES2"   NUTS("ES21","ES22","ES23","ES24","ES51","ES52","ES53","ES62")
    "ES3"   NUTS("ES30","ES41","ES42")
    "ES4"   NUTS("ES61","ES43")
    "FR1"   NUTS("FR81","FR82","FR83")
    "FR2"   NUTS("FR61","FR62","FR63")
    "FR3"   NUTS("FR51","FR52","FR53")
    "FR4"   NUTS("FR41","FR42","FR43")
    "FR5"   NUTS("FR10","FR21","FR22","FR23","FR24","FR25","FR26","FR30","FR71","FR72")
    "IT1"   NUTS("ITC1","ITC2","ITC3","ITC4","ITD1","ITD2","ITD3","ITD4","ITD5")
    "IT2"   NUTS("ITE1","ITE2","ITE3","ITE4")
    "IT3"   NUTS("ITF1","ITF2","ITF3","ITF4","ITF5","ITF6")
    "NO1"   NUTS("NO01","NO02","NO03","NO04","NO05")
    "NO2"   NUTS("NO06")
    "NO3"   NUTS("NO07")
    "PO1"   NUTS("PL21","PL22","PL51","PL52")
    "PO2"   NUTS("PL11","PL12","PL31","PL32","PL33","PL34")
    "PO3"   NUTS("PL41","PL42","PL43","PL61","PL62","PL63")
    "SE1"   NUTS("SE22")
    "SE2"   NUTS("SE22")
    "SE3"   NUTS("SE11","SE12","SE21","SE23","SE31")
    "SE4"   NUTS("SE33")
    "UK1"   NUTS("UKC1","UKC2","UKD1","UKD2","UKD3","UKD4","UKD5","UKE1","UKE2","UKE3","UKE4","UKF1","UKF2","UKF3","UKG1","UKG2","UKG3","UKH1","UKH2","UKH3","UKI1","UKI2","UKJ1","UKJ2","UKJ3","UKJ4","UKL1","UKL2")
    "UK2"   NUTS("UKM2","UKM3","UKM5","UKM6")
    "UK3"   NUTS("UKN0")
]

# Test region to see if mixed GADM/NUTS notation works
const testreg = [
    "NO"     GADM("Norway")
    "DK"     NUTS("DK")
    "ISL"    (NUTS("SE214"), GADM("Åland"), GADM(["Sweden","Kalmar"], "Borgholm","Mörbylånga"))  # Gotland, Åland, Öland
]