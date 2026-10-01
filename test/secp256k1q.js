const fieldTests = require("./tester/fieldtests.js");

const bigInt = require("big-integer");

fieldTests(new bigInt("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16), "secp256k1q");
