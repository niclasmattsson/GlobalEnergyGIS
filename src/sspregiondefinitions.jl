# Aggregation on the five region level

const ssp5 = Dict(
	# R5.2OECD = Includes the OECD 90 and EU member states and candidates.
	"OECD" => [
		"Albania", "Australia", "Austria", "Belgium", "Bosnia and Herzegovina", "Bulgaria", "Canada",
		"Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany",
		"Greece", "Guam", "Hungary", "Iceland", "Ireland", "Italy", "Japan", "Latvia", "Lithuania",
		"Luxembourg", "Malta", "Montenegro", "Netherlands", "New Zealand", "Norway", "Poland",
		"Portugal", "Puerto Rico", "Romania", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
		"Switzerland", "Macedonia", "Turkey", "United Kingdom", "United States"
	],

	# R5.2REF = Countries from the Reforming Economies of Eastern Europe and the Former Soviet Union.
	"REF" => [
		"Armenia", "Azerbaijan", "Belarus", "Georgia", "Kazakhstan", "Kyrgyzstan", "Moldova",
		"Russia", "Tajikistan", "Turkmenistan", "Ukraine", "Uzbekistan"
	],

	# R5.2ASIA = The region includes most Asian countries with the exception of the Middle East, Japan
	# and Former Soviet Union states. China includes Hong Kong and Macao, excludes Taiwan.
	"ASIA" => [
		"Afghanistan", "Bangladesh", "Bhutan", "Brunei", "Cambodia", "China", "North Korea", "Fiji",
		"French Polynesia", "India", "Indonesia", "Laos", "Malaysia", "Maldives", "Micronesia",
		"Mongolia", "Myanmar", "Nepal", "New Caledonia", "Pakistan", "Papua New Guinea", "Philippines",
		"South Korea", "Samoa", "Singapore", "Solomon Islands", "Sri Lanka", "Taiwan", "Thailand",
		"Timor-Leste", "Vanuatu", "Vietnam"
	],

	# R5.2MAF = This region includes the countries of the Middle East and Africa.
	"MAF" => [
		"Algeria", "Angola", "Bahrain", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon",
		"Cape Verde", "Central African Republic", "Chad", "Comoros", "Republic of Congo", "Côte d'Ivoire",
		"Democratic Republic of the Congo", "Djibouti", "Egypt", "Equatorial Guinea", "Eritrea", "Ethiopia",
		"Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Iran", "Iraq", "Israel", "Jordan", "Kenya",
		"Kuwait", "Lebanon", "Lesotho", "Liberia", "Libya", "Madagascar", "Malawi", "Mali", "Mauritania",
		"Mauritius", "Mayotte", "Morocco", "Mozambique", "Namibia", "Niger", "Nigeria", "Palestina", "Oman",
		"Qatar", "Rwanda", "Reunion", "Saudi Arabia", "Senegal", "Sierra Leone", "Somalia", "South Africa",
		"South Sudan", "Sudan", "Swaziland", "Syria", "Togo", "Tunisia", "Uganda", "United Arab Emirates",
		"Tanzania", "Western Sahara", "Yemen", "Zambia", "Zimbabwe"
	],

	# R5.2LAM = This region includes the countries of Latin America and the Caribbean.
	"LAM" => [
		"Argentina", "Aruba", "Bahamas", "Barbados", "Belize", "Bolivia", "Brazil", "Chile", "Colombia",
		"Costa Rica", "Cuba", "Dominican Republic", "Ecuador", "El Salvador", "French Guiana", "Grenada",
		"Guadeloupe", "Guatemala", "Guyana", "Haiti", "Honduras", "Jamaica", "Martinique", "Mexico",
		"Nicaragua", "Panama", "Paraguay", "Peru", "Suriname", "Trinidad and Tobago",
		"Virgin Islands, U.S.", "Uruguay", "Venezuela"
	]
)



# Aggregation on the 32 region level

const ssp32 = Dict(
	"ANUZ" => ["Australia", "New Zealand"],
	"BRA" => ["Brazil"],
	"CAN" => ["Canada"],
	# Central Asia
	"CAS" => ["Armenia", "Azerbaijan", "Georgia", "Kazakhstan", "Kyrgyzstan", "Tajikistan",
				"Turkmenistan", "Uzbekistan"],
	# China (Mainland, Hongkong, Macao; excl. Taiwan)
	"CHN" => ["China", "Hong Kong", "Macao"],
	# Eastern Europe (excl. former Soviet Union and EU member states)
	"EEU" => ["Albania", "Bosnia and Herzegovina", "Croatia", "Montenegro", "Serbia", "Macedonia"],
	# Eastern Europe, former Soviet Union (excl. Russia and EU members)
	"EEU-FSU" => ["Belarus", "Moldova", "Ukraine"],
	"EFTA" => ["Iceland", "Norway", "Switzerland"],
	# New EU member states that joined as of 2004 - high income.
	"EU12-H" => ["Cyprus", "Czech Republic", "Estonia", "Hungary", "Malta", "Poland", "Slovakia", "Slovenia"],
	# New EU member states that joined as of 2004 - medium income.
	"EU12-M" => ["Bulgaria", "Latvia", "Lithuania", "Romania"],
	# This region includes European Union member states that joined prior to 2004.
	"EU15" => ["Austria", "Belgium", "Denmark", "Finland", "France", "Germany", "Greece", "Ireland", "Italy",
				"Luxembourg", "Netherlands", "Portugal", "Spain", "Sweden", "United Kingdom"],
	"IDN" => ["Indonesia"],
	"IND" => ["India"],
	"JPN" => ["Japan"],
	"KOR" => ["South Korea"],
	# This region includes the countries of Latin America (excl. Brazil, Mexico) - low income.
	"LAM-L" => ["Belize", "Guatemala", "Haiti", "Honduras", "Nicaragua"],
	# This region includes the countries of Latin America (excl. Brazil, Mexico) - medium and high income.
	# [Netherlands Antilles not in GADM (new status), replaced with islands listed after Venezuela]
	"LAM-M" => ["Antigua and Barbuda", "Argentina", "Bahamas", "Barbados", "Bermuda", "Bolivia", "Chile",
				"Colombia", "Costa Rica", "Cuba", "Dominica", "Dominican Republic", "Ecuador", "El Salvador",
				"French Guiana", "Grenada", "Guadeloupe", "Guyana", "Jamaica", "Martinique",
				"Panama", "Paraguay", "Peru", "Saint Kitts and Nevis", "Saint Lucia",
				"Saint Vincent and the Grenadines", "Suriname", "Trinidad and Tobago", "Uruguay", "Venezuela",
				"Aruba", "Bonaire, Sint Eustatius and Saba", "Sint Maarten", "Curaçao"],
	# This region includes the countries of Middle East Asia - high income.
	"MEA-H" => ["Bahrain", "Israel", "Kuwait", "Oman", "Qatar", "Saudi Arabia", "United Arab Emirates"],
	# This region includes the countries of Middle East Asia - low and medium income.
	"MEA-M" => ["Iran", "Iraq", "Jordan", "Lebanon", "Palestina", "Syria", "Yemen"], 
	"MEX" => ["Mexico"],
	# This region includes the countries of North Africa.
	"NAF" => ["Algeria", "Egypt", "Libya", "Morocco", "Tunisia", "Western Sahara"],
	# This region includes the countries of Other Asia - former Centrally Planned Asia.
	"OAS-CPA" => ["Cambodia", "Laos", "Mongolia", "Vietnam"],
	# This region includes the countries of Other Asia - low income.
	"OAS-L" => ["Bangladesh", "North Korea", "Fiji", "Micronesia", "Myanmar", "Nepal", "Papua New Guinea",
				"Philippines", "Samoa", "Solomon Islands", "Timor-Leste", "Tonga", "Vanuatu"],
	# This region includes the countries of Other Asia - medium and high income.
	"OAS-M" => ["Bhutan", "Brunei", "French Polynesia", "Guam", "Malaysia", "Maldives", "New Caledonia",
				"Singapore", "Sri Lanka", "Thailand"],
	"PAK" => ["Pakistan", "Afghanistan"],
	"RUS" => ["Russia"],
	"SAF" => ["South Africa"],
	# This region includes the countries of Subsahara Africa (excl. South Africa) - low income.
	"SSA-L" => ["Benin", "Burkina Faso", "Burundi", "Cameroon", "Cape Verde", "Central African Republic",
				"Chad", "Comoros", "Republic of Congo", "Côte d'Ivoire", "Democratic Republic of the Congo",
				"Djibouti", "Eritrea", "Ethiopia", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya",
				"Lesotho", "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania", "Mozambique", "Niger",
				"Nigeria", "Rwanda", "São Tomé and Príncipe", "Senegal", "Sierra Leone", "Somalia",
				"South Sudan", "Sudan", "Swaziland", "Togo", "Uganda", "Tanzania", "Zambia", "Zimbabwe"],
	# This region includes the countries of Subsahara Africa (excl. South Africa) - medium and high income.
	"SSA-M" => ["Angola", "Botswana", "Equatorial Guinea", "Gabon", "Mauritius", "Mayotte", "Namibia",
				"Reunion", "Seychelles"],
	"TUR" => ["Turkey"],
	"TWN" => ["Taiwan"],
	"USA" => ["United States", "Puerto Rico", "Virgin Islands, U.S."]
)