

$(document).keyup(function(event) {
    if ($("#textInput.geneName").is(":focus") && (event.key == "Enter")) {
        $("#Button.geneName").click();
    }
});


