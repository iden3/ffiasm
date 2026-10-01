const fieldTests = require("./tester/fieldtests.js");

const bigInt = require("big-integer");

fieldTests(new bigInt("21888242871839275222246405745257275088696311157297823662689037894645226208583"), "bn128q");
