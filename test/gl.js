const fieldTests = require("./tester/fieldtests.js");

const bigInt = require("big-integer");

fieldTests(new bigInt("FFFFFFFF00000001", 16), "gl");
