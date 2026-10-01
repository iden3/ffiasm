const fieldTests = require("./tester/fieldtests.js");

const bigInt = require("big-integer");

fieldTests(new bigInt("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16), "secp256k1r");
