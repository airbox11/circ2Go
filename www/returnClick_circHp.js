$(document).keyup(function(event) {
    if ($("#textInput.geneName.heatmap").is(":focus") && (event.keyCode == 13)) {
        $("#Button.geneName.heatmap").click();
    }
});
