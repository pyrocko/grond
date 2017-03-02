function changeRundir(rundir) {
        var form = document.createElement("form");

        //Move the submit function to another variable
        //so that it doesn't get overwritten.
        form._submit_function_ = form.submit;

        form.setAttribute("method", "post");
        form.setAttribute("action", window.location.href);

        var rundirField = document.createElement("input");
        rundirField.setAttribute("type", "hidden");
        rundirField.setAttribute("name", "rundir");
        rundirField.setAttribute("value", rundir);

        form.appendChild(rundirField);

        document.body.appendChild(form);
        form._submit_function_(); //Call the renamed function.
    }
