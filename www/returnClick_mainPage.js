$(document).keyup(function(event) {
    if ($("#textInput.geneName").is(":focus") && (event.keyCode == 13)) {
        $("#mpActB1").click();
    }
});
