#' @title exportAgilentIL
#' @description Export RHermes inclusion lists as Agilent
#' Auto MS/MS Preferred/Exclude tables.
#'
#' @param struct RHermesExp object.
#' @param id Inclusion List ID.
#' @param file Output basename.
#' @param maxOver Maximum number of overlapping IL entries.
#' @param sepFiles Whether to export one file per injection.
#' @param maxInjections Maximum number of injections.
#' @param ppm Precursor mass tolerance (ppm).
#' @param charge Charge state.
#' @param isolationWidth One of "Narrow", "Medium", or "Wide".
#'
#' @export

setGeneric(
    "exportAgilentIL",
    function(struct,
             id,
             file = "./AgilentPreferred",
             maxOver = 5,
             sepFiles = TRUE,
             maxInjections = 9999,
             ppm = 20,
             charge = 1,
             isolationWidth = "Narrow") {

        standardGeneric("exportAgilentIL")
    }
)

setMethod(
    "exportAgilentIL",
    c("RHermesExp", "numeric"),

function(struct,
         id,
         file = "./AgilentPreferred",
         maxOver = 5,
         sepFiles = TRUE,
         maxInjections = 9999,
         ppm = 20,
         charge = 1,
         isolationWidth = "Narrow") {

    validObject(struct)

    if (length(struct@data@MS2Exp) == 0) {
        stop("This object doesn't have any ILs")
    }

    if (!between(id, 1, length(struct@data@MS2Exp))) {
        stop("Please enter a valid IL number")
    }

    isolationWidth <- match.arg(
        isolationWidth,
        c("Narrow", "Medium", "Wide")
    )

    isolationWidth <- switch(
        isolationWidth,
        Narrow = "Narrow (~1.3 m/z)",
        Medium = "Medium (~4 m/z)",
        Wide = "Wide (~8 m/z)"
    )

    IL <- struct@data@MS2Exp[[id]]@IL@IL

    plan <- injectionPlanner(
        IL,
        maxOver = maxOver,
        byMaxInt = TRUE,
        injections = maxInjections,
        mode = "continuous"
    )

    convert_plan <- function(p) {

        rt_center <- ((p$start + p$end) / 2) / 60
        rt_window <- ((p$end - p$start) / 2) / 60

        data.frame(
            On = "True",
            `Prec. m/z` = round(p$mass, 5),
            `Delta m/z (ppm)` = ppm,
            Z = charge,
            `Prec. Type` = "Preferred",
            `Ret. Time (min)` = round(rt_center, 2),
            `Delta Ret. Time (min)` = round(rt_window, 2),
            `Iso. Width` = isolationWidth,
            `Collision Energy` = "",
            check.names = FALSE
        )
    }

    write_agilent <- function(df, outfile) {

        con <- file(outfile, open = "wt")

        writeLines(
            "AutoPreferredExcludeMSMSTable",
            con
        )

        write.table(
            df,
            file = con,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )

        close(con)
    }

    if (sepFiles) {

        for (i in seq_along(plan)) {

            out <- convert_plan(plan[[i]])

            write_agilent(
                out,
                paste0(file, "_Injection_", i, ".csv")
            )
        }

    } else {

        out <- do.call(
            rbind,
            lapply(plan, convert_plan)
        )

        write_agilent(
            out,
            paste0(file, "_complete.csv")
        )
    }

    invisible(NULL)
})